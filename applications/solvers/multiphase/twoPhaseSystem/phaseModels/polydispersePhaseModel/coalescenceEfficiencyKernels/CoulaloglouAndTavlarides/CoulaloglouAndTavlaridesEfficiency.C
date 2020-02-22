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

#include "CoulaloglouAndTavlaridesEfficiency.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace coalescenceEfficiencyKernels
{
    defineTypeNameAndDebug(CoulaloglouAndTavlarides, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyKernel,
        CoulaloglouAndTavlarides,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::
CoulaloglouAndTavlarides
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& continuousPhase
)
:
    coalescenceEfficiencyKernel(dict, mesh, continuousPhase),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    Ceff_("Ceff", dimless, dict),
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
    ),
    muf_
    (
        IOobject
        (
            "CoulaloglouAndTavlarides:muf",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", dimDynamicViscosity, 0.0)
    )
{
    Ceff_.dimensions().reset(inv(sqr(dimLength)));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::
~CoulaloglouAndTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::update
(
    const fluidThermo& thermo,
    const turbulenceModel& turb
)
{
    epsilonf_ = turb.epsilon();
    epsilonf_.max(SMALL);
    muf_ = thermo.mu();
}


Foam::scalar Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::Pc
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

    scalar rho = fluid_.phase2().rho()[celli];
    scalar sigma = fluid_.sigma().value();

    return
        Foam::exp
        (
          - Ceff_.value()
           *muf_[celli]*epsilonf_[celli]*rho/sqr(sigma)
           *pow4(cbrt(v1)*cbrt(v2)/(cbrt(v1) + cbrt(v2)))
        );
}

// ************************************************************************* //
