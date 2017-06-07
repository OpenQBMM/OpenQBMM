/*---------------------------------------------------------------------------*\
 *  =========                 |
 *  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 *   \\    /   O peration     |
 *    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
 *     \\/     M anipUgation  |
 * -------------------------------------------------------------------------------
 * License
 *    This file is part of OpenFOAM.
 *
 *    OpenFOAM is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General PUplic License as pUplished by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    OpenFOAM is distributed in the hope that it will be usefUg, but WITHOUT
 *    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *    FITNESS FOR A PARTICUgAR PURPOSE.  See the GNU General PUplic License
 *    for more details.
 *
 *    You shoUgd have received a copy of the GNU General PUplic License
 *    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \*---------------------------------------------------------------------------*/


#include "anisotropicGaussianModel.H"
#include "fvm.H"
#include "fvc.H"
#include "fixedValueFvPatchFields.H"
#include "wallFvPatch.H"
#include "emptyFvPatch.H"
#include "coupledFvPatch.H"
#include "twoPhaseSystem.H"


// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

void Foam::RASModels::anisotropicGaussianModel::solveDilute
(
    const surfaceScalarField& h2f
)
{

	const dimensionedScalar& deltaT = mesh_.time().deltaT();

    // Lookup non-const refrences to alpha and U so they can be updated
    volScalarField& alphaRef =
        mesh_.lookupObjectRef<volScalarField>
        (
            IOobject::groupName("alpha", phase_.name())
        );

    volVectorField& URef =
        mesh_.lookupObjectRef<volVectorField>
        (
            IOobject::groupName("U", phase_.name())
        );

	tmp<surfaceScalarField> th1f(1.0 - h2f);
	const surfaceScalarField& h1f = th1f();

	Pp_ = Theta_*symmTensor::I - Sigma_;
	Pp_.correctBoundaryConditions();

    // Set old moment values
	volScalarField M0_old(phase_);
	volVectorField M1_old(phase_*phase_.U());
	volSymmTensorField M2_old(phase_*(Pp_ + sqr(phase_.U())));

    // Predictor step
	calcMomentFluxes(h1f);

	volScalarField M0("M0", M0_old - 0.5*fvc::surfaceIntegrate(F0_)*deltaT);
	volVectorField M1("M1", M1_old - 0.5*fvc::surfaceIntegrate(F1_)*deltaT);
	volSymmTensorField M2("M2", M2_old - 0.5*fvc::surfaceIntegrate(F2_)*deltaT);

	M0.max(SMALL);

    // Update alpha, U, and particle pressure
    alphaRef = M0;
    alphaRef.correctBoundaryConditions();

    URef = M1/M0;
    URef.correctBoundaryConditions();

    Pp_ = M2/M0 - sqr(phase_.U());
    Pp_.correctBoundaryConditions();

    // Corrector step
	calcMomentFluxes(h1f);

	M0 = M0_old - fvc::surfaceIntegrate(F0_)*deltaT;
	M1 = M1_old - fvc::surfaceIntegrate(F1_)*deltaT;
	M2 = M2_old - fvc::surfaceIntegrate(F2_)*deltaT;

	M0.max(SMALL);

	alphaRef = M0;
    alphaRef.correctBoundaryConditions();

	URef = M1/M0;
    URef.correctBoundaryConditions();

	Pp_ = M2/M0 - sqr(phase_.U());

	forAll(Pp_,i)
	{
		if(Pp_[i].xx() < SMALL) Pp_[i].xx() = SMALL;
		if(Pp_[i].yy() < SMALL) Pp_[i].yy() = SMALL;
		if(Pp_[i].zz() < SMALL) Pp_[i].zz() = SMALL;
	}
    Pp_.correctBoundaryConditions();

	Theta_ = 1.0/3.0*tr(Pp_);
	Sigma_ = Theta_*symmTensor::I - Pp_;

	Theta_.max(0);
	Theta_.min(100);

	Theta_.correctBoundaryConditions();
	Sigma_.correctBoundaryConditions();
}


void Foam::RASModels::anisotropicGaussianModel::calcMomentFluxes
(
    const surfaceScalarField& h1f
)
{

	const fvPatchList& patches = mesh_.boundary();
	const surfaceVectorField& Sf = mesh_.Sf();

	F0_ == dimensionedScalar("zero", F0_.dimensions(), 0.0);
	F1_ == dimensionedVector("zero", F1_.dimensions(), vector::zero);
	F2_ == dimensionedSymmTensor("zero", F2_.dimensions(), symmTensor::zero);

	if(cellInv_)
	{

		const labelUList& owner = mesh_.owner();

		forAll(phase_, celli)
		{
			const cell& cProp = mesh_.cells()[celli];

			quadrature_.calcHermiteQuadrature(phase_.U()[celli],Pp_[celli]);

			scalarField hwl(phase_[celli]*weights_);

			// Loop over all the current cell faces to calculate the flux
			forAll(cProp,locFacei)
			{
				const  label& facei = cProp[locFacei];

				if(mesh_.isInternalFace(facei))
				{

					scalarField facePhi(abscissae_ & Sf[facei]);
					// Sf is always pointing away from owner cell
					if (celli == owner[facei])
						facePhi *= pos(facePhi);
					else
						facePhi *= neg(facePhi);

					scalarField wPhi(hwl*facePhi);

					F0_[facei] += sum(wPhi);
					F1_[facei] += sum(wPhi*abscissae_);
					F2_[facei] += sum(wPhi*sqr(abscissae_));

				}
				else
				{
					register label patchi =
                        mesh_.boundaryMesh().whichPatch(facei);
					const fvPatch& currPatch = patches[patchi];

					if ( ! isA<emptyFvPatch> ( currPatch ) )
					{
						register label pFacei = facei - currPatch.start();
						const vector& sf = currPatch.Sf()[pFacei];

						//  outgoing flux on the wall
						scalarField facePhi( sf & abscissae_);
						facePhi *= pos(facePhi);
						scalarField wPhi(hwl*facePhi);

						scalar sFlux = sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] =
                            sum(wPhi*abscissae_);
						F2_.boundaryFieldRef()[patchi][pFacei] =
                            sum(wPhi*sqr(abscissae_));

						// refletion incoming flux on the wall
						if ( isA<wallFvPatch> ( patches[patchi] ) )
						{

							F0_.boundaryFieldRef()[patchi][pFacei] = 0.0;

							//  reflection flux on the wall
							vector nf(sf/mag(sf));
							vectorField Uw
							(
                                abscissae_
                              - (1.0 + ew_)*(abscissae_ & nf)*nf
                            );
							wPhi *= -(1.0 - phiw_);

							F1_.boundaryFieldRef()[patchi][pFacei] +=
                                sum(wPhi*Uw);
							symmTensor tFlux(sum(wPhi*sqr(Uw)));
							F2_.boundaryFieldRef()[patchi][pFacei] += tFlux;

							scalar dfls_mf = phiw_*sFlux;
							scalar dfls_ef = phiw_*tr(tFlux);

							if(dfls_mf*dfls_ef > SMALL )
							{

								scalarField facePhio
								(
                                    sf
                                  & quadrature_.hermiteOriginalAbscissas()
                                );
								facePhio *= neg(facePhio);
								scalarField wPhio(hwl*facePhio);

								scalar sFluxo = mag(sum(wPhio));
								scalar tFluxo = mag
                                (
                                    sum
                                    (
                                        wPhio
                                       *magSqr
                                        (
                                            quadrature_.hermiteOriginalAbscissas()
                                        )
                                    )
                                );

								scalar sig =
                                    sqrt(dfls_ef*sFluxo/(dfls_mf*tFluxo));
								scalar pp = dfls_mf/sFluxo/sig;

								vectorField  Ud
								(
                                    sig*quadrature_.hermiteOriginalAbscissas()
                                );
								scalarField  facePhid( sf & Ud );
								facePhid *= neg(facePhid);
								scalarField  wPhid(pp*hwl*facePhid);

								F1_.boundaryFieldRef()[patchi][pFacei] +=
                                    sum(wPhid*Ud);
								F2_.boundaryFieldRef()[patchi][pFacei] +=
                                    sum(wPhid*sqr(Ud));
							}

						}

					}
				}

			}

		}

		forAll(patches, patchi)
		{
			const fvPatch& currPatch = patches[patchi];

			if
            (
               !(
                   isA<wallFvPatch> ( currPatch )
                || isA<emptyFvPatch> ( currPatch )
                )
            )
			{

				const vectorField& Sfbf = currPatch.Sf();

				// incoming flux from the periodic boundary
				if ( isA<coupledFvPatch> (currPatch) )
				{

					scalarField alphaNf
					(
                        phase_.boundaryField()[patchi].patchNeighbourField()
                    );
					vectorField UpNf
					(
                        phase_.U().boundaryField()[patchi].patchNeighbourField()
                    );
					symmTensorField PpNf
					(
                        Pp_.boundaryField()[patchi].patchNeighbourField()
                    );

					forAll(Sfbf,pFacei)
					{
						quadrature_.calcHermiteQuadrature
						(
						    UpNf[pFacei],
						    PpNf[pFacei]
						);

						scalarField facePhi(Sfbf[pFacei] & abscissae_);
						facePhi *= neg(facePhi);
						scalarField wPhi(alphaNf[pFacei]*weights_*facePhi);

						F0_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] +=
                            sum(wPhi*abscissae_);
						F2_.boundaryFieldRef()[patchi][pFacei] +=
                            sum(wPhi*sqr(abscissae_));

					}

				}
				else
				{
					// incoming flux from the inlet boundary
					const scalarField& alphaPf(phase_.boundaryField()[patchi]);
					const vectorField& UpPf(phase_.U().boundaryField()[patchi]);
					const symmTensorField& PpPf(Pp_.boundaryField()[patchi]);

					forAll(Sfbf,pFacei)
					{

						quadrature_.calcHermiteQuadrature
						(
						    UpPf[pFacei],
						    PpPf[pFacei]
						);

						scalarField facePhi(Sfbf[pFacei] & abscissae_);
						facePhi *= neg(facePhi);
						scalarField wPhi(alphaPf[pFacei]*weights_*facePhi);

						F0_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] +=
                            sum(wPhi*abscissae_);
						F2_.boundaryFieldRef()[patchi][pFacei] +=
                            sum(wPhi*sqr(abscissae_));

					}
				}
			}
		}

		F0_ *= h1f;
		F1_ *= h1f;
		F2_ *= h1f;

	}
	else
	{

		forAll(ownnei_,di)
		{

			surfaceScalarField alphaf
			(
                fvc::interpolate(phase_,ownnei_[di], "alphap")
            );
			surfaceVectorField Upf(fvc::interpolate(phase_.U(),ownnei_[di],"Up"));
			surfaceSymmTensorField Ppf(fvc::interpolate(Pp_,ownnei_[di],"Pp"));

			alphaf *= h1f;

			forAll(F0_, facei)
			{



				quadrature_.calcHermiteQuadrature(Upf[facei], Ppf[facei]);

				scalarField facePhi(abscissae_ & Sf[facei]);

				// Sf is always pointing away from the owner cell to neightbor cell
				if(di == 0 )
					facePhi *= pos(facePhi);
				else
					facePhi *= neg(facePhi);

				scalarField wPhi(alphaf[facei]*weights_*facePhi);

				scalar sWphi(sum(wPhi));
				if(sWphi != 0.0)
				{
					F0_[facei] += sWphi;
					F1_[facei] += sum(wPhi*abscissae_);
					F2_[facei] += sum(wPhi*sqr(abscissae_));
				}

			}


			forAll(patches, patchi)
			{

				const fvPatch& currPatch = patches[patchi];
				const vectorField& Sfbf = currPatch.Sf();

				const scalarField& alphaPf(alphaf.boundaryField()[patchi]);
				const vectorField& UpPf(Upf.boundaryField()[patchi]);
				const symmTensorField& PpPf(Ppf.boundaryField()[patchi]);

				// flux from the periodic boundary
				if ( isA<coupledFvPatch> (currPatch) )
				{

					forAll(Sfbf,pFacei)
					{
						quadrature_.calcHermiteQuadrature
						(
						    UpPf[pFacei],
						    PpPf[pFacei]
						);

						scalarField facePhi(Sfbf[pFacei] & abscissae_);
						if(di == 0 )
							facePhi *= pos(facePhi);
						else
							facePhi *= neg(facePhi);

						scalarField wPhi(alphaPf[pFacei]*weights_*facePhi);

						scalar sWphi(sum(wPhi));
						if(sWphi != 0.0)
						{
							F0_.boundaryFieldRef()[patchi][pFacei] += sWphi;
							F1_.boundaryFieldRef()[patchi][pFacei] +=
                                sum(wPhi*abscissae_);
							F2_.boundaryFieldRef()[patchi][pFacei] +=
                                sum(wPhi*sqr(abscissae_));
						}

					}

				}
				else if(isA<wallFvPatch> ( currPatch ) && (di == 0))
				{
					// wall boundary only need to solve once

					F0_.boundaryFieldRef()[patchi] == 0.0;

					forAll(Sfbf,pFacei)
					{
						const vector& sf = Sfbf[pFacei];

						scalarList hwl(alphaPf[pFacei]*weights_);

						quadrature_.calcHermiteQuadrature
						(
						    UpPf[pFacei],
						    PpPf[pFacei]
						);

						//  outgoing flux on the wall
						scalarField facePhi( sf & abscissae_);
						facePhi *= pos(facePhi);
						scalarField wPhi(hwl*facePhi);

						scalar sFlux = sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] =
                            sum(wPhi*abscissae_);
						F2_.boundaryFieldRef()[patchi][pFacei] =
                            sum(wPhi*sqr(abscissae_));

						//  reflection flux on the wall
						vector nf(sf/mag(sf));
						vectorField Uw
						(
                            abscissae_
                          - (1.0 + ew_)*(abscissae_ & nf)*nf
                        );
						wPhi *= -(1.0 - phiw_);

						F1_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*Uw);
						symmTensor tFlux(sum(wPhi*sqr(Uw)));
						F2_.boundaryFieldRef()[patchi][pFacei] += tFlux;

						scalar dfls_mf = phiw_*sFlux;
						scalar dfls_ef = phiw_*tr(tFlux);

						if(dfls_mf*dfls_ef > SMALL )
						{

							scalarField facePhio
							(
                                sf & quadrature_.hermiteOriginalAbscissas()
                            );

							facePhio *= neg(facePhio);
							scalarField wPhio(hwl*facePhio);

							scalar sFluxo = mag(sum(wPhio));
							scalar tFluxo = mag
							(
                                sum
                                (
                                    wPhio
                                   *magSqr
                                    (
                                        quadrature_.hermiteOriginalAbscissas()
                                    )
                                )
                            );


							scalar sig = sqrt(dfls_ef*sFluxo/(dfls_mf*tFluxo));
							scalar pp = dfls_mf/sFluxo/sig;

							vectorField  Ud
							(
                                sig*quadrature_.hermiteOriginalAbscissas()
                            );

							scalarField  facePhid( sf & Ud );
							facePhid *= neg(facePhid);
							scalarField  wPhid(pp*hwl*facePhid);

							F1_.boundaryFieldRef()[patchi][pFacei] +=
                                sum(wPhid*Ud);
							F2_.boundaryFieldRef()[patchi][pFacei] +=
                                sum(wPhid*sqr(Ud));
						}

					}

				}
				else if((!isA<emptyFvPatch>(currPatch)) && (di == 0) )
				{
					//  inlet/outlet boundary only need to solve once
					//  fixed value boundary condition, i.e. inlet only consider
                    //  incoming flux otherwise, zeroGradient outlet type, only
                    //  consider outgoing flux

					bool inlet =
                        isA<fixedValueFvPatchScalarField>
                        (
                            phase_.boundaryField()[patchi]
                        );

					forAll(Sfbf,pFacei)
					{


						quadrature_.calcHermiteQuadrature
						(
						    UpPf[pFacei],
						    PpPf[pFacei]
						);

						scalarField facePhi(Sfbf[pFacei] & abscissae_);

						if(inlet)
							facePhi *= neg(facePhi);
						else
							facePhi *= pos(facePhi);

						scalarField wPhi(alphaPf[pFacei]*weights_*facePhi);

						F0_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi*abscissae_);
						F2_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi*sqr(abscissae_));

					}
				}
			}
		}
	}
	return ;
}


void Foam::RASModels::anisotropicGaussianModel::updateh2Fn()
{

    const dimensionedScalar smallPpk
    (
        "small",
        dimensionSet(0, 2, -2, 0, 0),
        SMALL
    );

    gs0_ = radialModel_->g0(phase_, alphaMinFriction_, alphaMax_);

    // This calculates the h2 function for the dense regime transport.
    if(h2FnMethod_.match("alphaG0"))
    {
        h2Fn_ = 1.0 - 1.0/(1.0 + sqr(phase_)*pow(gs0_,h2FnParaPow_));
    }
    else
    {
        if(h2FnMethod_.match("particlePressure"))
        {
            volScalarField ppk(max(phase_*Theta_, smallPpk));
            volScalarField pps(4.0*eta_*phase_*gs0_*ppk + ppfr_);
            h2Fn_ = pow(pps/(pps + ppk), h2FnParaPow_);
        }
        else
        {
            FatalErrorIn("kineticTheoryModel::updateh2Fn: invalid h2FnMethod") << abort(FatalError);
        }

    }
    h2Fn_.correctBoundaryConditions();
}


Foam::scalar Foam::RASModels::anisotropicGaussianModel::maxUxDx()
{
    scalar maxUxDx = 0.0;
    tmp<volScalarField>  tVv;

	tVv =
        (mag(phase_.U().component(vector::X))
      + quadrature_.maxAbs()*sqrt(Pp_.component(symmTensor::XX)))
	   /pDxyz_.component(vector::X);
	maxUxDx = gMax(tVv());
	tVv.clear();

	tVv =
        (mag(phase_.U().component(vector::Y))
      + quadrature_.maxAbs()*sqrt(Pp_.component(symmTensor::YY)))
       /pDxyz_.component(vector::Y);
	maxUxDx = max(maxUxDx, gMax(tVv()));
	tVv.clear();

	tVv =
        (mag(phase_.U().component(vector::Z))
      + quadrature_.maxAbs()*sqrt(Pp_.component(symmTensor::ZZ)))
	   /pDxyz_.component(vector::Z);
	maxUxDx = max(maxUxDx, gMax(tVv()));
	tVv.clear();

    return maxUxDx;

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::anisotropicGaussianModel::anisotropicGaussianModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    kineticTheoryModel
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName,
        type
    ),
    mesh_(alpha.mesh()),
	cellInv_(coeffDict_.lookupOrDefault("cellInv", true)),
	eta_(0.5*(1.0 + e_)),
	Sigma_
    (
       IOobject
       (
           IOobject::groupName("Sigma", phase_.name()),
           mesh_.time().timeName(),
           mesh_,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       devRhoReff()/phase_.rho()
    ),
	Pp_
	(
	    IOobject
	    (
           IOobject::groupName("Pp", phase_.name()),
	        mesh_.time().timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    Theta_*symmTensor::I - Sigma_,
	    Sigma_.boundaryField().types()
	),
    ppfr_
    (
       IOobject
       (
           IOobject::groupName("ppfr", phase_.name()),
           mesh_.time().timeName(),
           mesh_,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh_,
       dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
    ),
    h2Fn_
    (
        IOobject
        (
            IOobject::groupName("h2Fn", phase_.name()),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        1.0
    ),
    h2FnMethod_(coeffDict_.lookup("h2FnMethod")),
    h2FnParaPow_(readLabel(coeffDict_.lookup("h2FnParaPow"))),
	quadrature_(3, coeffDict_.lookupOrDefault("nNodePerDim", 4)),
	weights_(quadrature_.hermiteWeights()),
	abscissae_(quadrature_.hermiteAbscissas()),
	ew_(coeffDict_.lookupOrDefault("wallRestitutionCoefficient", 1.0)),
	phiw_(coeffDict_.lookupOrDefault("wallSpecularityCoefficient", 0.0)),
	ownnei_(2),
	F0_
	(
	    IOobject
	    (
	        "F0",
	        mesh_.time().timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimensionedScalar("zero", dimless*dimVol/dimTime, scalar(0))
	),
	F1_
	(
	    IOobject
	    (
	        "F1",
	        mesh_.time().timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimensionedVector
	    (
            "zero",
            dimVelocity*dimVol/dimTime,
            Zero
        )
	),
	F2_
	(
	    IOobject
	    (
	        "dM2",
	        mesh_.time().timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimensionedSymmTensor
	    (
            "zero",
            pow(dimVelocity,2)*dimVol/dimTime,
            Zero
        )
	),
    pDxyz_
	(
		IOobject
		(
			"dxyz",
			mesh_.time().timeName(),
			mesh_
		),
		mesh_,
		dimensionedVector("zero", dimless, Zero)
	)
{
	if (!cellInv_)
	{
		forAll(ownnei_,i)
		{
			ownnei_.set
			(
			    i,
			    new surfaceScalarField
			    (
			        IOobject
			        (
			            "phi",
			            mesh_.time().timeName(),
			            mesh_,
			            IOobject::NO_READ,
			            IOobject::NO_WRITE
			        ),
			        mesh_,
			        0.0
			    )
			);
		}

		ownnei_[0] == 1.0;
		ownnei_[1] == -1.0;
	}

	const faceList & ff = mesh_.faces();
	const pointField & pp = mesh_.points();

	forAll (mesh_.C(), celli)
	{
		const cell & cc = mesh_.cells()[celli];
		labelList pLabels(cc.labels(ff));
		pointField pLocal(pLabels.size(), vector::zero);

		forAll(pLabels, pointi) pLocal[pointi] = pp[pLabels[pointi]];

		pDxyz_[celli][0] =
            Foam::max(pLocal & vector(1,0,0))
          - Foam::min(pLocal & vector(1,0,0));
		pDxyz_[celli][1] =
            Foam::max(pLocal & vector(0,1,0))
          - Foam::min(pLocal & vector(0,1,0));
		pDxyz_[celli][2] =
            Foam::max(pLocal & vector(0,0,1))
          - Foam::min(pLocal & vector(0,0,1));
	}

	forAll(mesh_.boundary(), patchi)
    {
        pDxyz_.boundaryFieldRef()[patchi] == vector(1,1,1);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::anisotropicGaussianModel::~anisotropicGaussianModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::RASModels::anisotropicGaussianModel::pPrime() const
{

    return
        (
            h2Fn_
          + 4.0*phase_*eta_
           *(
                2.0*radialModel_->g0(phase_, alphaMinFriction_, alphaMax_)
              + radialModel_->g0prime
                (
                    volScalarField(phase_),
                    alphaMinFriction_,
                    alphaMax_
                )*phase_
            )
        )*Theta_
      + frictionalStressModel_->frictionalPressurePrime
        (
            phase_,
            alphaMinFriction_,
            alphaMax_
        );

}


void Foam::RASModels::anisotropicGaussianModel::correct()
{
    if (dilute_)
    {
        updateh2Fn();
        surfaceScalarField h2Fnf = fvc::interpolate(h2Fn_);
        solveDilute(h2Fnf);
        surfaceScalarField& phi =
            mesh_.lookupObjectRef<surfaceScalarField>
            (
                IOobject::groupName
                (
                    "phi",
                    phase_.name()
                )
            );
        phi = fvc::flux(phase_.U());
        dilute_ = false;
    }
    else
    {
        // Local references
        const dimensionedScalar smallRT("small",dimensionSet(0, 0, -1, 0, 0), SMALL);
        const volScalarField& alpha = phase_;
        const volScalarField& rho = phase_.rho();
        const surfaceScalarField& alphaPhi = phase_.alphaPhi();
        const volScalarField& d = phase_.d();
        const scalar sqrtPi = sqrt(constant::mathematical::pi);

        // Particle Strain-rate tensor
        tmp<volTensorField> tgradU(fvc::grad(phase_.U()));
        const volTensorField& gradU(tgradU());
        volSymmTensorField D(symm(gradU));  // D = 0.5*(gradU + gradU^T)
        volSymmTensorField Sp(D - (1.0/3.0)*tr(D)*I);

        // particle continuity error
        volScalarField contErr =
            fvc::ddt(alpha)
          + fvc::div(alphaPhi);

        // Radial distribution model
        gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

        volScalarField alphaSqr(sqr(alpha));

        volScalarField thetaSqrt(sqrt(Theta_));

        // bulk viscosity
        lambda_ = (8.0/3.0)/sqrtPi*d*eta_*alphaSqr*gs0_*thetaSqrt;

        // particle pressure, move bulk viscosity part to tau
        volScalarField  PsCoeff = alpha*(h2Fn_ + 4.0*eta_*alpha*gs0_);

        volScalarField rTaupAlpha
        (
            "rTaupAlpha",
            phase_.fluid().Kd()/rho
          + smallRT
        );

        volScalarField rTaucAlphaCoeff((6.0/sqrtPi/d)*gs0_*max(alphaSqr, residualAlpha_));
        volScalarField rTaucAlpha("rTaucAlpha", rTaucAlphaCoeff*thetaSqrt);

        // Particle viscosity
        volScalarField nutCoeff
        (
            0.5*alphaSqr*(1.0 + (8.0/5.0)*eta_*(3*eta_ - 2.0)*alpha*gs0_)
           *h2Fn_*( 1.0 + (8.0/5.0)*eta_*alpha*gs0_)
        );

        nut_ =
            nutCoeff*Theta_/(rTaupAlpha + eta_*(2.0-eta_)*rTaucAlpha)
          + (3.0/5.0)*lambda_;

        volSymmTensorField tau(2.0*nut_*Sp + lambda_*tr(D)*I);

        // 'thermal' conductivity
        volScalarField kappaCoeff
        (
            2.5*alphaSqr*(1.0 + (12.0/5.0)*sqr(eta_)
           *(4.0*eta_ - 3.0)*alpha*gs0_)
           *h2Fn_*(1.0 + (12.0/5.0)*eta_*alpha*gs0_)
        );
        kappa_ =
            kappaCoeff*Theta_
           /(3.0*rTaupAlpha + 4.0*eta_*(41.0-33.0*eta_)*rTaucAlpha)
          + (3.0/2.0)*lambda_;

        // AG corrections
        // update collisional pressure
        volScalarField ppc(4.0*eta_*alphaSqr*gs0_*Theta_ - lambda_*tr(D)) ;

        volSymmTensorField S2flux(2.0*((h2Fn_*alpha*Theta_ + ppc)*Sp - nut_*( twoSymm(Sp & gradU) - (2.0/3.0)*(Sp && gradU)*I)));

        fvSymmTensorMatrix SigmaEqn
        (
            fvm::ddt(alpha, Sigma_)
          + fvm::div(alphaPhi, Sigma_)
          - fvc::Sp(contErr, Sigma_)

          - fvm::laplacian(2.0/3.0*kappa_, Sigma_)
         ==
            S2flux
          + fvm::Sp
            (
                -(2.0*rTaupAlpha + (3.0-e_)*(1.0+e_)/2.0*rTaucAlpha),
                Sigma_
            )

        );

        SigmaEqn.relax();
        SigmaEqn.solve();

        Sigma_.correctBoundaryConditions();

        // Construct the granular temperature equation
        fvScalarMatrix ThetaEqn
        (

            fvm::ddt(alpha, Theta_)
          + fvm::div(alphaPhi, Theta_)
          - fvc::Sp(contErr, Theta_)
          - fvm::laplacian(2.0/3.0*kappa_, Theta_)
         ==
            fvm::SuSp(- 2.0/3.0*((PsCoeff*I) && gradU), Theta_)
          + 2.0/3.0*( tau && gradU)
          + fvm::Sp(- (2.0*rTaupAlpha + (1.0 - sqr(e_))*rTaucAlpha), Theta_)
        );

        ThetaEqn.relax();
        ThetaEqn.solve();

        Theta_.max(0);
        Theta_.min(100);

        Theta_.correctBoundaryConditions();

        thetaSqrt = sqrt(Theta_);

        // update bulk viscosity
        lambda_ = (8.0/3.0)/sqrtPi*d*eta_*alphaSqr*gs0_*thetaSqrt;

        // update Particle viscosity
        nut_ =
            nutCoeff*Theta_
           /(rTaupAlpha + eta_*(2.0-eta_)*rTaucAlphaCoeff*thetaSqrt)
          + (3.0/5.0)*lambda_;

        // Frictional pressure
        ppfr_ = frictionalStressModel_->frictionalPressure
        (
            phase_,
            alphaMinFriction_,
            alphaMax_
        );

        // Add frictional shear viscosity
        nut_ += frictionalStressModel_->nu
        (
            phase_,
            alphaMinFriction_,
            alphaMax_,
            ppfr_,
            Sp
        );

        // Limit viscosity
        nut_.min(maxNut_);

        if (debug)
        {
            Info<< "    max(Theta) = " << max(Theta_).value() << nl
                << "    min(Theta) = " << min(Theta_).value() << nl
                << "    max(nut) = " << max(nut_).value() << endl;
        }

        dilute_ = true;
    }
}
// ************************************************************************* //
