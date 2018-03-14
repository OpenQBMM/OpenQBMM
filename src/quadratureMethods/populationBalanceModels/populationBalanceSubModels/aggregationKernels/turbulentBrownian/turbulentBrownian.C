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

#include "turbulentBrownian.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(turbulentBrownian, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        turbulentBrownian,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::turbulentBrownian
::turbulentBrownian
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh),
    flThermo_(mesh_.lookupObject<fluidThermo>(basicThermo::dictName)),
    flTurb_
    (
        mesh_.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    ),
    T_(flThermo_.T()),
    rho_(flThermo_.rho()),
    mu_(flThermo_.mu()),
    epsilon_(flTurb_.epsilon())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::turbulentBrownian
::~turbulentBrownian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::turbulentBrownian::Ka
(
    const scalar& abscissa1,
    const scalar& abscissa2,
    const label celli
) const
{
    return 2.0*Foam::constant::physicoChemical::k.value()*T_[celli]
        *sqr(abscissa1 + abscissa2)/(3.0*mu_[celli]
        *max(abscissa1*abscissa2, SMALL)) + 4.0/3.0*pow3(abscissa1 + abscissa2)
        *sqrt(3.0*Foam::constant::mathematical::pi*epsilon_[1]
        /(10.0*mu_[celli]/rho_[celli]));
}

// ************************************************************************* //
