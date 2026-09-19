/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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
    const surfaceScalarField& phi,
    const supportType& support
)
:
    PDFTransportModel(name, dict, mesh),
    quadrature_(name, mesh, List<supportType>(1, support)),
    momentAdvection_
    (
        univariateMomentAdvection::New
        (
            quadrature_.subDict("momentAdvection"),
            quadrature_,
            phi,
            support
        )
    ),
    SSPRK2_(false)
{
    checkRealizableTransport(quadrature_.moments());

    // The explicit Euler step keeps the moments realizable under the Courant
    // limit of the advection scheme but is of first order in time, which
    // takes the accuracy of a second-order scheme such as zeta away at the
    // Courant numbers it allows. The second-order strong stability
    // preserving Runge-Kutta scheme is two such steps and their average
    // with the moments of the start of the step: a convex combination of
    // realizable moment sets, so it is realizable under the same limit.
    const word timeIntegration
    (
        quadrature_.subDict("momentAdvection").lookupOrDefault<word>
        (
            "timeIntegration",
            "Euler"
        )
    );

    if (timeIntegration == "SSPRK2")
    {
        SSPRK2_ = true;
    }
    else if (timeIntegration != "Euler")
    {
        FatalIOErrorInFunction(quadrature_.subDict("momentAdvection"))
            << "Unknown timeIntegration " << timeIntegration << nl
            << "    Valid ones are Euler and SSPRK2." << nl
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFTransportModel
::~univariatePDFTransportModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::univariatePDFTransportModel::solveSSPRK2()
{
    volScalarMomentFieldSet& moments = quadrature_.moments();
    const scalar deltaT = mesh_.time().deltaTValue();

    // The moments at the start of the step, which every stage, and every
    // outer iteration of a step, starts from
    PtrList<scalarField> startMoments(moments.size());

    forAll(moments, mi)
    {
        startMoments.set
        (
            mi,
            new scalarField(moments[mi].oldTime().primitiveField())
        );
    }

    // First stage: an Euler step from the start of the step
    momentAdvection_().update();

    forAll(moments, mi)
    {
        volScalarMoment& m = moments[mi];

        m.primitiveFieldRef() =
            startMoments[mi]
          - deltaT*momentAdvection_().divMoments()[mi].primitiveField();

        m.correctBoundaryConditions();
    }

    // The advection schemes build the moments of the faces from the moments
    // of the cells, so the quadrature is not inverted between the stages

    // Second stage: an Euler step from the first, averaged with the start
    momentAdvection_().update();

    forAll(moments, mi)
    {
        volScalarMoment& m = moments[mi];

        const scalarField secondStage
        (
            m.primitiveField()
          - deltaT*momentAdvection_().divMoments()[mi].primitiveField()
        );

        m.primitiveFieldRef() = 0.5*(startMoments[mi] + secondStage);

        m.correctBoundaryConditions();
    }

    // The implicit source, the diffusion of the moments, from the advected
    // moments: the time derivative the matrix carries is taken from the
    // start of the step, so the advection over the step is added to it as a
    // source, and the equation says (m - m_advected)/deltaT = source(m)
    forAll(moments, mi)
    {
        volScalarMoment& m = moments[mi];

        volScalarField advection
        (
            IOobject
            (
                "advection",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(m.dimensions()/dimTime, Zero)
        );

        advection.primitiveFieldRef() =
            (m.primitiveField() - startMoments[mi])/deltaT;

        fvScalarMatrix momentEqn
        (
            fvm::ddt(m)
          - advection
         ==
            implicitMomentSource(m)
        );

        momentEqn.relax();
        momentEqn.solve();
    }
}


void Foam::PDFTransportModels::univariatePDFTransportModel::solve()
{
    if (SSPRK2_)
    {
        solveSSPRK2();

        quadrature_.updateQuadrature();

        if (solveMomentSources())
        {
            this->explicitMomentSource();
        }

        return;
    }

    momentAdvection_().update();

    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(quadrature_.nMoments());

    // Solve moment transport equations
    forAll(quadrature_.moments(), momenti)
    {
        volScalarMoment& m = quadrature_.moments()[momenti];

        momentEqns.set
        (
            momenti,
            new fvScalarMatrix
            (
                fvm::ddt(m)
              + momentAdvection_().divMoments()[momenti]
              ==
                implicitMomentSource(m)
            )
        );
    }

    forAll (momentEqns, mEqni)
    {
        momentEqns[mEqni].relax();
        momentEqns[mEqni].solve();
    }

    quadrature_.updateQuadrature();

    if (solveMomentSources())
    {
        this->explicitMomentSource();
    }
}


// ************************************************************************* //
