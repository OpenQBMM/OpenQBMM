/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 Alberto Passalacqua
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

#include "CoulaloglouAndTavlaridesFrequency.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace coalescenceFrequencyKernels
{
   defineTypeNameAndDebug(CoulaloglouAndTavlarides, 0);

    addToRunTimeSelectionTable
    (
        coalescenceFrequencyKernel,
        CoulaloglouAndTavlarides,
        dictionary
    );

} // End namespace coalescenceFrequencyKernels
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::
CoulaloglouAndTavlarides
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& continuousPhase
)
:
    coalescenceFrequencyKernel(dict, mesh, continuousPhase),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    epsilonf_
    (
        IOobject
        (
            "CoulaloglouAndTavlarides:epsilonf",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", sqr(dimVelocity)/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::
~CoulaloglouAndTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::aggregationKernels::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::update
(
    const fluidThermo& thermo,
    const turbulenceModel& turb
)
{
    epsilonf_ = turb.epsilon();
    epsilonf_.max(SMALL);
}


Foam::scalar Foam::populationBalanceSubModels::aggregationKernels::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::omega
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli
) const
{
    scalar pi = constant::mathematical::pi;
    scalar v1 = pow3(d1)*pi/6.0;
    scalar v2 = pow3(d2)*pi/6.0;

    return
        (pow(v1, 2.0/3.0) + pow(v2, 2.0/3.0))
       *sqrt(pow(v1, 2.0/9.0) + pow(v2, 2.0/9.0))
       *cbrt(epsilonf_[celli])/(1.0 + fluid_.phase1()[celli]);
}

// ************************************************************************* //
