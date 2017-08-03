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
    name_(name),
    mixingKernel_
    (
        Foam::mixingSubModels::mixingKernel::New
        (
            dict.subDict("mixingKernel")
        )
    ),
    diffusionModel_
    (
        Foam::mixingSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::mixingModels::turbulentMixing::~turbulentMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<fvScalarMatrix> Foam::PDFTransportModels::mixingModels
::turbulentMixing::momentDiffusion
(
    const volUnivariateMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::mixingModels::turbulentMixing
::phaseSpaceConvection
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> gSource
    (
        new volScalarField
        (
            IOobject
            (
                "gSource",
                phi_.mesh().time().timeName(),
                phi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phi_.mesh(),
            dimensionedScalar("zero", moment.dimensions()/dimTime, 0.0)
        )
    );

    return gSource;
}

void Foam::PDFTransportModels::mixingModels::turbulentMixing
::explicitMomentSource()
{}

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::mixingModels::turbulentMixing::implicitMomentSource
(
    const volUnivariateMoment& moment
)
{
    const volUnivariateMomentFieldSet& moments = (*this).quadrature().moments();

    return mixingKernel_->K(moment, moments);
}

Foam::scalar
Foam::PDFTransportModels::mixingModels::turbulentMixing::cellMomentSource
(
    label& momentOrder,
    label& celli
)
{
    return 0.0;
}

Foam::scalar Foam::PDFTransportModels::mixingModels::turbulentMixing::realizableCo
()
{
    return univariatePDFTransportModel::realizableCo();
}

void Foam::PDFTransportModels::mixingModels::turbulentMixing::solve
()
{
    univariatePDFTransportModel::solve();
}

// ************************************************************************* //
