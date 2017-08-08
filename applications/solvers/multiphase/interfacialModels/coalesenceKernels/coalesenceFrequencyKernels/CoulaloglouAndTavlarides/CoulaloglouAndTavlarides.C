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
namespace coalesenceFrequencyKernels
{
    defineTypeNameAndDebug(CoulaloglouAndTavlarides, 0);

    addToRunTimeSelectionTable
    (
        coalesenceFrequencyKernel,
        CoulaloglouAndTavlarides,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::coalesenceFrequencyKernels::
CoulaloglouAndTavlarides::CoulaloglouAndTavlarides
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalesenceFrequencyKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    epsilon_(fluid_.phase2().turbulence().epsilon()())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::coalesenceFrequencyKernels::
CoulaloglouAndTavlarides::~CoulaloglouAndTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::coalesenceFrequencyKernels::
CoulaloglouAndTavlarides::omega
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    return
        cbrt(epsilon_)*sqr(d1 + d2)
       *sqrt(pow(d1, 2.0/3.0) + pow(d2, 2.0/3.0));
}

// ************************************************************************* //
