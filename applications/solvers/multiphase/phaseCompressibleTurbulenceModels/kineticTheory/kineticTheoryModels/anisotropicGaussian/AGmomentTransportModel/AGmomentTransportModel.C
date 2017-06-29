/*---------------------------------------------------------------------------*\
 *  =========                 |
 *  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 *   \\    /   O peration     |
 *    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
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


#include "AGmomentTransportModel.H"
#include "fvc.H"
#include "fixedValueFvPatchFields.H"
#include "wallFvPatch.H"
#include "emptyFvPatch.H"
#include "coupledFvPatch.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AGmomentTransportModel::AGmomentTransportModel
(
    const fvMesh& mesh,
    const dictionary& dict,
	const phaseModel& phase,
    volScalarField& Theta,
    volSymmTensorField& Sigma
)
:
	cellInv_(dict.lookupOrDefault("cellInv", true)),
	mesh_(mesh),
	phase_(phase),
	alphap_
	(
        mesh.lookupObjectRef<volScalarField>
        (
            IOobject::groupName("alpha", phase.name())
        )
    ),
	Up_
	(
        mesh.lookupObjectRef<volVectorField>
        (
            IOobject::groupName("U", phase.name())
        )
    ),
	Theta_(Theta),
	Sigma_(Sigma),
	Pp_
	(
	    IOobject
	    (
            IOobject::groupName("Pp", phase.name()),
	        mesh_.time().timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
	    Theta_*symmTensor::I - Sigma_,
	    Sigma_.boundaryField().types()
	),
	hq_(3, dict.lookupOrDefault("nHerNodePerDim", 4)),
	hqWeigs_(hq_.hermiteWeights()),
	hqAbsc_(hq_.hermiteAbscissas()),
	ew_(dict.lookupOrDefault("wallRestitutionCoefficient", 1.0)),
	phiw_(dict.lookupOrDefault("wallSpecularityCoefficient", 0.0)),
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
	    dimensionedVector("zero", dimVelocity*dimVol/dimTime, vector::zero)
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
	    dimensionedSymmTensor("zero", pow(dimVelocity,2)*dimVol/dimTime, symmTensor::zero)
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
		vector::zero
	)
{
	if( !cellInv_)
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

	const faceList & ff = mesh.faces();
	const pointField & pp = mesh.points();

	forAll ( mesh.C(), celli)
	{
		const cell & cc = mesh_.cells()[celli];
		labelList pLabels(cc.labels(ff));
		pointField pLocal(pLabels.size(), vector::zero);

		forAll(pLabels, pointi) pLocal[pointi] = pp[pLabels[pointi]];

		pDxyz_[celli][0] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
		pDxyz_[celli][1] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
		pDxyz_[celli][2] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
	}

	forAll(mesh_.boundary(), patchi) pDxyz_.boundaryFieldRef()[patchi] == vector(1,1,1);


}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::AGmomentTransportModel::~AGmomentTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar  Foam::AGmomentTransportModel::maxUxDx() const
{
    scalar maxUxDx = 0.0;
    tmp<volScalarField>  tVv;

	tVv = (mag(Up_.component(vector::X)) + hq_.maxAbs()*sqrt(Pp_.component(symmTensor::XX)))
		  /pDxyz_.component(vector::X);
	maxUxDx = gMax(tVv());
	tVv.clear();

	tVv = (mag(Up_.component(vector::Y)) + hq_.maxAbs()*sqrt(Pp_.component(symmTensor::YY)))
		  /pDxyz_.component(vector::Y);
	maxUxDx = max(maxUxDx, gMax(tVv()));
	tVv.clear();

	tVv = (mag(Up_.component(vector::Z)) + hq_.maxAbs()*sqrt(Pp_.component(symmTensor::ZZ)))
		  /pDxyz_.component(vector::Z);
	maxUxDx = max(maxUxDx, gMax(tVv()));
	tVv.clear();

    return maxUxDx;
}


void Foam::AGmomentTransportModel::solve
(
    const surfaceScalarField& h2f
)
{

	const dimensionedScalar& deltaT = mesh_.time().deltaT();

	tmp<surfaceScalarField> th1f(1.0 - h2f);
	const surfaceScalarField& h1f = th1f();

	Pp_ = Theta_*symmTensor::I - Sigma_;
	Pp_.correctBoundaryConditions();

	volScalarField M0_old(alphap_);
	volVectorField M1_old(alphap_*Up_);
	volSymmTensorField M2_old(alphap_*(Pp_ + sqr(Up_)));

	calcMomentFluxes(h1f);

	volScalarField M0("M0", M0_old - 0.5*fvc::surfaceIntegrate(F0_)*deltaT);
	volVectorField M1("M1", M1_old - 0.5*fvc::surfaceIntegrate(F1_)*deltaT);
	volSymmTensorField M2("M2", M2_old - 0.5*fvc::surfaceIntegrate(F2_)*deltaT);

	M0.max(SMALL);
	alphap_ = M0;
	Up_ = M1/M0;
	Pp_ = M2/M0 - sqr(Up_);

	alphap_.correctBoundaryConditions();
	Up_.correctBoundaryConditions();
	Pp_.correctBoundaryConditions();

	calcMomentFluxes(h1f);

	M0 = M0_old - fvc::surfaceIntegrate(F0_)*deltaT;
	M1 = M1_old - fvc::surfaceIntegrate(F1_)*deltaT;
	M2 = M2_old - fvc::surfaceIntegrate(F2_)*deltaT;

	M0.max(SMALL);
	alphap_ = M0;
	Up_ = M1/M0;
	Pp_ = M2/M0 - sqr(Up_);

	forAll(Pp_,i)
	{
		if(Pp_[i].xx() < SMALL) Pp_[i].xx() = SMALL;
		if(Pp_[i].yy() < SMALL) Pp_[i].yy() = SMALL;
		if(Pp_[i].zz() < SMALL) Pp_[i].zz() = SMALL;
	}

	alphap_.correctBoundaryConditions();
	Up_.correctBoundaryConditions();
	Pp_.correctBoundaryConditions();

	Theta_ = 1.0/3.0*tr(Pp_);
	Sigma_ = Theta_*symmTensor::I - Pp_;

	Theta_.max(0);
	Theta_.min(100);

	Theta_.correctBoundaryConditions();
	Sigma_.correctBoundaryConditions();

	surfaceScalarField& phip =
        mesh_.lookupObjectRef<surfaceScalarField>
        (
            IOobject::groupName("phi", phase_.name())
        );
    phip = fvc::flux(Up_);
}


void Foam::AGmomentTransportModel::calcMomentFluxes
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

		forAll(alphap_, celli)
		{
			const cell& cProp = mesh_.cells()[celli];

			hq_.calcHermiteQuadrature(Up_[celli],Pp_[celli]);

			scalarField hwl(alphap_[celli]*hqWeigs_);

			// Loop over all the current cell faces to calculate the flux
			forAll(cProp,locFacei)
			{
				const  label& facei = cProp[locFacei];

				if(mesh_.isInternalFace(facei))
				{

					scalarField facePhi(hqAbsc_ & Sf[facei]);
					// Sf is always pointing away from owner cell
					if(celli == owner[facei])
						facePhi *= pos(facePhi);
					else
						facePhi *= neg(facePhi);

					scalarField wPhi(hwl*facePhi);

					F0_[facei] += sum(wPhi);
					F1_[facei] += sum(wPhi*hqAbsc_);
					F2_[facei] += sum(wPhi*sqr(hqAbsc_));

				}
				else
				{
					register label patchi = mesh_.boundaryMesh().whichPatch(facei);
					const fvPatch& currPatch = patches[patchi];

					if ( ! isA<emptyFvPatch> ( currPatch ) )
					{
						register label pFacei = facei - currPatch.start();
						const vector& sf = currPatch.Sf()[pFacei];

						//  outgoing flux on the wall
						scalarField facePhi( sf & hqAbsc_);
						facePhi *= pos(facePhi);
						scalarField wPhi(hwl*facePhi);

						scalar sFlux = sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi*hqAbsc_);
						F2_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi*sqr(hqAbsc_));

						// refletion incoming flux on the wall
						if ( isA<wallFvPatch> ( patches[patchi] ) )
						{

							F0_.boundaryFieldRef()[patchi][pFacei] = 0.0;

							//  reflection flux on the wall
							vector nf(sf/mag(sf));
							vectorField Uw(hqAbsc_ - (1.0 + ew_)*(hqAbsc_ & nf)*nf);
							wPhi *= -(1.0 - phiw_);

							F1_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*Uw);
							symmTensor tFlux(sum(wPhi*sqr(Uw)));
							F2_.boundaryFieldRef()[patchi][pFacei] += tFlux;

							scalar dfls_mf = phiw_*sFlux;
							scalar dfls_ef = phiw_*tr(tFlux);

							if(dfls_mf*dfls_ef > SMALL )
							{

								scalarField facePhio(sf & hq_.hermiteOriginalAbscissas());
								facePhio *= neg(facePhio);
								scalarField wPhio(hwl*facePhio);

								scalar sFluxo = mag(sum(wPhio));
								scalar tFluxo = mag(sum(wPhio*magSqr(hq_.hermiteOriginalAbscissas())));

								scalar sig = std::sqrt(dfls_ef*sFluxo/(dfls_mf*tFluxo));
								scalar pp = dfls_mf/sFluxo/sig;

								vectorField  Ud(sig*hq_.hermiteOriginalAbscissas());
								scalarField  facePhid( sf & Ud );
								facePhid *= neg(facePhid);
								scalarField  wPhid(pp*hwl*facePhid);

								F1_.boundaryFieldRef()[patchi][pFacei] += sum(wPhid*Ud);
								F2_.boundaryFieldRef()[patchi][pFacei] += sum(wPhid*sqr(Ud));
							}

						}

					}
				}

			}

		}

		forAll(patches, patchi)
		{
			const fvPatch& currPatch = patches[patchi];

			if (!(isA<wallFvPatch> ( currPatch ) || isA<emptyFvPatch> ( currPatch )) )
			{

				const vectorField& Sfbf = currPatch.Sf();

				// incoming flux from the periodic boundary
				if ( isA<coupledFvPatch> (currPatch) )
				{

					scalarField alphaNf(alphap_.boundaryField()[patchi].patchNeighbourField());
					vectorField UpNf(Up_.boundaryField()[patchi].patchNeighbourField());
					symmTensorField PpNf(Pp_.boundaryField()[patchi].patchNeighbourField());

					forAll(Sfbf,pFacei)
					{
						hq_.calcHermiteQuadrature
						(
						    UpNf[pFacei],
						    PpNf[pFacei]
						);

						scalarField facePhi(Sfbf[pFacei] & hqAbsc_);
						facePhi *= neg(facePhi);
						scalarField wPhi(alphaNf[pFacei]*hqWeigs_*facePhi);

						F0_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*hqAbsc_);
						F2_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*sqr(hqAbsc_));

					}

				}
				else
				{
					// incoming flux from the inlet boundary
					const scalarField& alphaPf(alphap_.boundaryField()[patchi]);
					const vectorField& UpPf(Up_.boundaryField()[patchi]);
					const symmTensorField& PpPf(Pp_.boundaryField()[patchi]);

					forAll(Sfbf,pFacei)
					{

						hq_.calcHermiteQuadrature
						(
						    UpPf[pFacei],
						    PpPf[pFacei]
						);

						scalarField facePhi(Sfbf[pFacei] & hqAbsc_);
						facePhi *= neg(facePhi);
						scalarField wPhi(alphaPf[pFacei]*hqWeigs_*facePhi);

						F0_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*hqAbsc_);
						F2_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*sqr(hqAbsc_));

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

			surfaceScalarField alphaf(fvc::interpolate(alphap_,ownnei_[di], "alphap"));
			surfaceVectorField Upf(fvc::interpolate(Up_,ownnei_[di],"Up"));
			surfaceSymmTensorField Ppf(fvc::interpolate(Pp_,ownnei_[di],"Pp"));

			alphaf *= h1f;

			forAll(F0_, facei)
			{



				hq_.calcHermiteQuadrature(Upf[facei], Ppf[facei]);

				scalarField facePhi(hqAbsc_ & Sf[facei]);

				// Sf is always pointing away from the owner cell to neightbor cell
				if(di == 0 )
					facePhi *= pos(facePhi);
				else
					facePhi *= neg(facePhi);

				scalarField wPhi(alphaf[facei]*hqWeigs_*facePhi);

				scalar sWphi(sum(wPhi));
				if(sWphi != 0.0)
				{
					F0_[facei] += sWphi;
					F1_[facei] += sum(wPhi*hqAbsc_);
					F2_[facei] += sum(wPhi*sqr(hqAbsc_));
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
						hq_.calcHermiteQuadrature
						(
						    UpPf[pFacei],
						    PpPf[pFacei]
						);

						scalarField facePhi(Sfbf[pFacei] & hqAbsc_);
						if(di == 0 )
							facePhi *= pos(facePhi);
						else
							facePhi *= neg(facePhi);

						scalarField wPhi(alphaPf[pFacei]*hqWeigs_*facePhi);

						scalar sWphi(sum(wPhi));
						if(sWphi != 0.0)
						{
							F0_.boundaryFieldRef()[patchi][pFacei] += sWphi;
							F1_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*hqAbsc_);
							F2_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*sqr(hqAbsc_));
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

						scalarList hwl(alphaPf[pFacei]*hqWeigs_);

						hq_.calcHermiteQuadrature
						(
						    UpPf[pFacei],
						    PpPf[pFacei]
						);

						//  outgoing flux on the wall
						scalarField facePhi( sf & hqAbsc_);
						facePhi *= pos(facePhi);
						scalarField wPhi(hwl*facePhi);

						scalar sFlux = sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi*hqAbsc_);
						F2_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi*sqr(hqAbsc_));

						//  reflection flux on the wall
						vector nf(sf/mag(sf));
						vectorField Uw(hqAbsc_ - (1.0 + ew_)*(hqAbsc_ & nf)*nf);
						wPhi *= -(1.0 - phiw_);

						F1_.boundaryFieldRef()[patchi][pFacei] += sum(wPhi*Uw);
						symmTensor tFlux(sum(wPhi*sqr(Uw)));
						F2_.boundaryFieldRef()[patchi][pFacei] += tFlux;

						scalar dfls_mf = phiw_*sFlux;
						scalar dfls_ef = phiw_*tr(tFlux);

						if(dfls_mf*dfls_ef > SMALL )
						{

							scalarField facePhio(sf & hq_.hermiteOriginalAbscissas());
							facePhio *= neg(facePhio);
							scalarField wPhio(hwl*facePhio);

							scalar sFluxo = mag(sum(wPhio));
							scalar tFluxo = mag(sum(wPhio*magSqr(hq_.hermiteOriginalAbscissas())));


							scalar sig = std::sqrt(dfls_ef*sFluxo/(dfls_mf*tFluxo));
							scalar pp = dfls_mf/sFluxo/sig;

							vectorField  Ud(sig*hq_.hermiteOriginalAbscissas());
							scalarField  facePhid( sf & Ud );
							facePhid *= neg(facePhid);
							scalarField  wPhid(pp*hwl*facePhid);

							F1_.boundaryFieldRef()[patchi][pFacei] += sum(wPhid*Ud);
							F2_.boundaryFieldRef()[patchi][pFacei] += sum(wPhid*sqr(Ud));
						}

					}

				}
				else if((!isA<emptyFvPatch>(currPatch)) && (di == 0) )
				{
					// inlet/outlet boundary only need to solve once
					// fixed value boundary condition, i.e. inlet only consider incoming flux
					// otherwise, zeroGradient outlet type,  only consider outgoing flux

					bool inlet = isA<fixedValueFvPatchScalarField>(alphap_.boundaryField()[patchi]);

					forAll(Sfbf,pFacei)
					{


						hq_.calcHermiteQuadrature
						(
						    UpPf[pFacei],
						    PpPf[pFacei]
						);

						scalarField facePhi(Sfbf[pFacei] & hqAbsc_);

						if(inlet)
							facePhi *= neg(facePhi);
						else
							facePhi *= pos(facePhi);

						scalarField wPhi(alphaPf[pFacei]*hqWeigs_*facePhi);

						F0_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi);
						F1_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi*hqAbsc_);
						F2_.boundaryFieldRef()[patchi][pFacei] = sum(wPhi*sqr(hqAbsc_));

					}
				}
			}
		}
	}

	return;
}


// ************************************************************************* //
