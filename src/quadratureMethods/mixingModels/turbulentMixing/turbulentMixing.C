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
    Copyright (C) 2019-2024 Alberto Passalacqua
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

#include "turbulentMixing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace mixingModels
{
    defineTypeNameAndDebug(turbulentMixing, 0);
    addToRunTimeSelectionTable
    (
        mixingModel,
        turbulentMixing,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::mixingModels::turbulentMixing::turbulentMixing
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    univariatePDFTransportModel(name, dict, phi.mesh(), phi, "01"),
    mixingModel(name, dict, phi),
    odeType(phi.mesh(), dict),
    name_(name),
    mixingKernel_
    (
        Foam::mixingSubModels::mixingKernel::New
        (
            dict.subDict("mixingKernel"),
            phi.mesh(),
            (*this).quadrature().moments()
        )
    ),
    diffusionModel_
    (
        Foam::mixingSubModels::mixingDiffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::mixingModels::turbulentMixing::~turbulentMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::PDFTransportModels::mixingModels::turbulentMixing
::explicitMomentSource()
{
    odeType::solve(quadrature_, 0);
}


bool Foam::PDFTransportModels::mixingModels::turbulentMixing
::solveMomentSources() const
{
    return odeType::solveSources_;
}


bool Foam::PDFTransportModels::mixingModels::turbulentMixing
::solveMomentOde() const
{
    return odeType::solveOde_;
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::mixingModels::turbulentMixing
::implicitMomentSource(const volScalarMoment& moment)
{
    return diffusionModel_->momentDiff(moment);
}


Foam::scalar Foam::PDFTransportModels::mixingModels::turbulentMixing
::realizableCo()
{
    return univariatePDFTransportModel::realizableCo();
}


void Foam::PDFTransportModels::mixingModels::turbulentMixing::solve()
{
    univariatePDFTransportModel::solve();
}


void
Foam::PDFTransportModels::mixingModels::turbulentMixing
::updateCellMomentSource(const label)
{}


Foam::scalar
Foam::PDFTransportModels::mixingModels::turbulentMixing
::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
)
{
    scalar source(0);

    source +=
        mixingKernel_->mixingSource
        (
            momentOrder[0],
            celli,
            0
        );

    return source;
}

// ************************************************************************* //
