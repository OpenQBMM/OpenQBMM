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

#include "oneQuarterMassRatio.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace daughterDistributions
{
    defineTypeNameAndDebug(oneQuarterMassRatio, 0);

    addToRunTimeSelectionTable
    (
        daughterDistribution,
        oneQuarterMassRatio,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::daughterDistributions::oneQuarterMassRatio
::oneQuarterMassRatio
(
    const dictionary& dict
)
:
    daughterDistribution(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::daughterDistributions::oneQuarterMassRatio
::~oneQuarterMassRatio()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::daughterDistributions::oneQuarterMassRatio::mD
(
    const label& order,
    const scalar& abscissa
) const
{
    scalar exponent = order/3.0;
    return
    (
        (pow(scalar(4), exponent) + 1.0)
      * pow(abscissa, order)
      / pow(scalar(5), exponent)
    );
}


Foam::scalar
Foam::populationBalanceSubModels::daughterDistributions::oneQuarterMassRatio
::mDMass
(
    const label& order,
    const scalar& abscissa
) const
{
    return
    (
        (pow(scalar(4), order) + 1.0)
      * pow(abscissa, order)
      / pow(label(5), order)
    );
}

// ************************************************************************* //
