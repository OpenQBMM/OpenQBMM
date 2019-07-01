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

#include "LehrEfficiency.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace coalescenceEfficiencyKernels
{
    defineTypeNameAndDebug(Lehr, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyKernel,
        Lehr,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Lehr::
Lehr
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& continuousPhase
)
:
    coalescenceEfficiencyKernel(dict, mesh, continuousPhase),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    sigma_(fluid_.sigma()),
    WeCrit_
    (
        dimensionedScalar::lookupOrDefault
        (
            "WeCrit",
            dict,
            dimVelocity,
            0.06
        )
    ),
    epsilonf_
    (
        IOobject
        (
            "Lehr:epsilonf",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", sqr(dimVelocity)/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Lehr::
~Lehr()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Lehr::update
(
    const fluidThermo& thermo,
    const turbulenceModel& turb
)
{
    epsilonf_ = turb.epsilon();
}


Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Lehr::Pc
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli
) const
{
    scalar dEq = 2.0/(1.0/d1 + 1.0/d2);

    return max
    (
        WeCrit_.value()*sigma_.value()/(fluid_.phase2().rho()[celli]*dEq)
       /max
        (
            sqrt(2.0)*pow(epsilonf_[celli]*sqrt(d1*d2), 1.0/3.0),
            mag(Ur)
        ),
        1.0
    );
}

// ************************************************************************* //
