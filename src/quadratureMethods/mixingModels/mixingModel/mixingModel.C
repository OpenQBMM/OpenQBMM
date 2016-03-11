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

#include "mixingModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::mixingModels::mixingModel::mixingModel
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    univariatePDFTransportModel(name, dict, U.mesh(), U, "01"),
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

Foam::PDFTransportModels::mixingModels::mixingModel::~mixingModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<fvScalarMatrix> Foam::PDFTransportModels::mixingModels
::mixingModel::momentDiffusion
(
    const volUnivariateMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::mixingModels::mixingModel
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
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    return gSource;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::mixingModels::mixingModel
::momentSource
(
    const volUnivariateMoment& moment
)
{
    const volUnivariateMomentFieldSet& moments = (*this).quadrature().moments();

    return mixingKernel_->K(moment, moments);
}

void Foam::PDFTransportModels::mixingModels
::mixingModel::solve
()
{
    univariatePDFTransportModel::solve();
}

// ************************************************************************* //
