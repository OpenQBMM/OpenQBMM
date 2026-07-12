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
    Copyright (C) 2019-2023 Alberto Passalacqua
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

#include "thompsonAggregation.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(thompsonAggregation, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        thompsonAggregation,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::thompsonAggregation
::thompsonAggregation
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::thompsonAggregation
::~thompsonAggregation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::thompsonAggregation::Ka
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli,
    const label environment
) const
{
    if
    (
        !std::isfinite(d1)
     || !std::isfinite(d2)
     || d1 < 0.0
     || d2 < 0.0
    )
    {
        return 0.0;
    }

    const scalar d13 = pow(d1, 3);
    const scalar d23 = pow(d2, 3);

    if (!std::isfinite(d13) || !std::isfinite(d23))
    {
        return 0.0;
    }

    const scalar ka =
        Ca_.value()
       *(
            pow(d13 - d23, 2)
           /(d13 + d23 + SMALL)
        );

    if (!std::isfinite(ka) || ka <= 0.0)
    {
        return 0.0;
    }

    return ka;
}

// ************************************************************************* //
