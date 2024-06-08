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
    k_
    (
        flTurb_.k()
    ),
    epsilon_
    (
        flTurb_.epsilon()
    ),
    moments_(moments)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernel::~mixingKernel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
