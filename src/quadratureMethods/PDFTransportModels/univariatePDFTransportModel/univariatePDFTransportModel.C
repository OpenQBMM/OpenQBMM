/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
    const surfaceScalarField& phi,
    const word& support
)
:
    PDFTransportModel(name, dict, mesh),
    name_(name),
    quadrature_(name, mesh, support),
    momentAdvection_
    (
        univariateMomentAdvection::New
        (
            quadrature_.subDict("momentAdvection"),
            quadrature_,
            phi,
            support
        )
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::univariatePDFTransportModel
::~univariatePDFTransportModel()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::univariatePDFTransportModel::solve()
{
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

    if (solveMomentSources())
    {
        this->explicitMomentSource();
    }

    forAll (momentEqns, mEqni)
    {
        volScalarField& m = quadrature_.moments()[mEqni];

        if (solveMomentSources())
        {
            momentEqns[mEqni] -= fvc::ddt(m);
        }

        momentEqns[mEqni].relax();
        momentEqns[mEqni].solve();
    }

    quadrature_.updateQuadrature();
}


// ************************************************************************* //
