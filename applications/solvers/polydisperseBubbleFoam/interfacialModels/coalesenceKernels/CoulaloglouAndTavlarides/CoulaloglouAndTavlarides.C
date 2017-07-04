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

#include "CoulaloglouAndTavlarides.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(CoulaloglouAndTavlarides, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        CoulaloglouAndTavlarides,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::CoulaloglouAndTavlarides
::CoulaloglouAndTavlarides
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    Ceff_(dict.lookup("Ceff")),
    sigma_(fluid_.sigma()),
    rho_(fluid_.phase2().rho()),
    mu_(fluid_.phase2().mu()),
    epsilon_(fluid_.phase2().turbulence().epsilon()())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::CoulaloglouAndTavlarides
::~CoulaloglouAndTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::CoulaloglouAndTavlarides::Ka
(
    const scalar& abscissa1,
    const scalar& abscissa2,
    const label& celli
) const
{
    return
        Ca_.value()*cbrt(epsilon_[celli])*sqr(abscissa1 + abscissa2)
       *sqrt(pow(abscissa1, 2.0/3.0) + pow(abscissa2, 2.0/3.0))
       *Foam::exp
        (
          - Ceff_.value()
           *sqrt
            (
                2.0*rho_[celli]*mu_[celli]*epsilon_[celli]/sqr(sigma_.value())
               *pow4(abscissa1*abscissa2/max(abscissa1 + abscissa2, SMALL))
            )
        );
}

// ************************************************************************* //
