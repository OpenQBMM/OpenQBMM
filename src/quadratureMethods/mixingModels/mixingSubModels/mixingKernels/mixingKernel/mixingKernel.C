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

#include "mixingKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
    defineTypeNameAndDebug(mixingKernel, 0);
    defineRunTimeSelectionTable(mixingKernel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernel::mixingKernel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volScalarMomentFieldSet& moments
)
:
    dict_(dict),
    mesh_(mesh),
    Cphi_
    (
        dict.lookupOrDefault
        (
            "Cphi",
            dimensionedScalar("CPhiDefault", dimless, 2.0)
        )
    ),
    Cmixing_
    (
        dict.lookupOrDefault
        (
            "Cmixing",
            dimensionedScalar("CmixingDefault", dimless, 1.0)
        )
    ),
    flTurb_
    (
        mesh_.lookupObject<turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                word::null
            )
        )
    ),
    k_(),
    epsilon_(),
    turbulenceTimeIndex_(-1),
    moments_(moments)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernel::~mixingKernel()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::mixingSubModels::mixingKernel::updateTurbulence() const
{
    // The kernels used to bind a reference to what the turbulence model
    // returned when they were built. A model that solves for k and epsilon
    // returns its own fields, which the reference followed; one that does
    // not, such as the k-omega family for epsilon, returns a field built for
    // the call, which the reference outlived.
    const label timeIndex = mesh_.time().timeIndex();

    if (timeIndex != turbulenceTimeIndex_)
    {
        k_ = flTurb_.k();
        epsilon_ = flTurb_.epsilon();
        turbulenceTimeIndex_ = timeIndex;
    }
}


const Foam::volScalarField& Foam::mixingSubModels::mixingKernel::k() const
{
    updateTurbulence();
    return k_();
}


const Foam::volScalarField&
Foam::mixingSubModels::mixingKernel::epsilon() const
{
    updateTurbulence();
    return epsilon_();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
