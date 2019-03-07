/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Alberto Passalacqua
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

#include "velocityPDFTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::velocityPDFTransportModel::velocityPDFTransportModel
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh,
    const word& support
)
:
    PDFTransportModel(name, dict, mesh),
    name_(name),
    quadrature_(name, mesh, support),
    momentAdvection_
    (
        velocityMomentAdvection::New
        (
            quadrature_.subDict("momentAdvection"),
            quadrature_,
            support
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::velocityPDFTransportModel::~velocityPDFTransportModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::velocityPDFTransportModel::solve()
{
    momentAdvection_().update();

    // List of moment transport equations
    PtrList<fvScalarMatrix> momentEqns(quadrature_.nMoments());

    // Solve moment transport equations
    forAll(quadrature_.moments(), momenti)
    {
        volVectorMoment& m = quadrature_.moments()[momenti];
        momentEqns.set
        (
            momenti,
            new fvScalarMatrix
            (
                fvm::ddt(m)
              + momentAdvection_().divMoments()[momenti]
            )
        );
    }

    if (solveMomentSources())
    {
        this->explicitMomentSource();
    }
    else
    {
        updateImplicitMomentSource();
    }

    forAll(quadrature_.moments(), mEqni)
    {
        volVectorMoment& m = quadrature_.moments()[mEqni];

        if (solveMomentSources())
        {
            momentEqns[mEqni] -= fvc::ddt(m);
        }
        else
        {
            // Solve moment transport excluding collisions
            momentEqns[mEqni].relax();
            momentEqns[mEqni].solve();

            //  Set moments.oldTime to moments transport is not neglected due to
            //  large collision source terms
            m.oldTime() = m;

            // Solve collisions
            momentEqns.set
            (
                mEqni,
                new fvScalarMatrix
                (
                    fvm::ddt(m)
                 ==
                    implicitMomentSource(m)
                )
            );
        }
    }

    forAll(quadrature_.moments(), mEqni)
    {
        momentEqns[mEqni].relax();
        momentEqns[mEqni].solve();
    }

    quadrature_.updateQuadrature();
}


// ************************************************************************* //
