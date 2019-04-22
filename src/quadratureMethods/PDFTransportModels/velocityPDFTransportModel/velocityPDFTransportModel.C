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

    // Solve moment transport equations
    forAll(quadrature_.moments(), momenti)
    {
        volVelocityMoment& m = quadrature_.moments()[momenti];
        fvScalarMatrix momentEqn
        (
            fvm::ddt(m)
          + momentAdvection_().divMoments()[momenti]
        );
        momentEqn.relax();
        momentEqn.solve();
    }
    quadrature_.updateQuadrature();

    if (solveMomentSources())
    {
        this->explicitMomentSource();
        updateImplicitMomentSource();


        bool update = false;
        forAll(quadrature_.moments(), mEqni)
        {
            const volVelocityMoment& m = quadrature_.moments()[mEqni];
            fvScalarMatrix iSource(implicitMomentSource(m));

            if (max(mag(iSource.source())) > small)
            {
                //  Set moments.oldTime so moments transport is not neglected due
                //  to large collision source terms
                quadrature_.moments()[mEqni].oldTime() = m;

                update = true;

                // Solve collisions
                fvScalarMatrix momentEqn
                (
                    fvm::ddt(m)
                 ==
                    implicitMomentSource(m)
                );

                momentEqn.relax();
                momentEqn.solve();
            }
        }

        if (update)
        {
            quadrature_.updateQuadrature();
        }
    }
}


// ************************************************************************* //
