/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
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
    updateImplicitMomentSource();
    
    forAll(quadrature_.moments(), momenti)
    {
        volVelocityMoment& m = quadrature_.moments()[momenti];
        fvScalarMatrix momentEqn
        (
            fvm::ddt(m)
          + momentAdvection_().divMoments()[momenti]
         ==
            implicitMomentSource(m)
        );
        momentEqn.relax();
        momentEqn.solve();
    }

    quadrature_.updateQuadrature();

    if (solveMomentSources())
    {
        this->explicitMomentSource();
    }
}


// ************************************************************************* //
