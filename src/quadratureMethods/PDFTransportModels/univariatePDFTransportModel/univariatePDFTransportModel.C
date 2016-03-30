/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "univariatePDFTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFTransportModel
::univariatePDFTransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const volVectorField& U,
    const word support
)
:
    PDFTransportModel(name, dict, mesh),
    name_(name),
    quadrature_(name, mesh, support),
    U_(U)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFTransportModel
::~univariatePDFTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::univariatePDFTransportModel
::updatePhysicalSpaceConvection
(
    surfaceScalarField& phiOwn,
    surfaceScalarField& phiNei
)
{
    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("nei", dimless, -1.0)
    );

    surfaceScalarField own
    (
        IOobject
        (
            "own",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("own", dimless, 1.0)
    );

    phiOwn = fvc::interpolate(U_, own, "reconstruct(U)") & mesh_.Sf();
    phiNei = fvc::interpolate(U_, nei, "reconstruct(U)") & mesh_.Sf();

    // Update interpolated nodes
    quadrature_.interpolateNodes();

    // Updated reconstructed moments
    quadrature_.momentsNei().update();
    quadrature_.momentsOwn().update();
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::univariatePDFTransportModel::physicalSpaceConvection
(
    const volUnivariateMoment& moment,
    const surfaceScalarField& phiOwn,
    const surfaceScalarField& phiNei
)
{
    dimensionedScalar zeroPhi("zero", phiNei.dimensions(), 0.0);

    tmp<volScalarField> divMoment
    (
        new volScalarField
        (
            IOobject
            (
                "divMoment",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    label order = moment.order();

    surfaceScalarField mFlux
    (
        quadrature_.momentsNei()[order]*min(phiNei, zeroPhi)
      + quadrature_.momentsOwn()[order]*max(phiOwn, zeroPhi)
    );

    fvc::surfaceIntegrate(divMoment.ref(), mFlux);
    divMoment.ref().dimensions().reset(moment.dimensions()/dimTime);

    return divMoment;
}

void Foam::PDFTransportModels::univariatePDFTransportModel::solve()
{
    surfaceScalarField phiOwn("phiOwn", fvc::interpolate(U_) & mesh_.Sf());
    surfaceScalarField phiNei("phiNei", phiOwn);
    updatePhysicalSpaceConvection(phiOwn, phiNei);

    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(quadrature_.nMoments());

    // Solve moment transport equations
    forAll(quadrature_.moments(), mI)
    {
        volUnivariateMoment& m = quadrature_.moments()[mI];

        momentEqns.set
        (
            mI,
            new fvScalarMatrix
            (
                fvm::ddt(m)
              + physicalSpaceConvection(m, phiOwn, phiNei)
              - momentDiffusion(m)
              ==
                momentSource(m)
              + phaseSpaceConvection(m)
            )
        );
    }

    forAll (momentEqns, mEqnI)
    {
        momentEqns[mEqnI].relax();
        momentEqns[mEqnI].solve();
    }

    quadrature_.updateQuadrature();
}


// ************************************************************************* //
