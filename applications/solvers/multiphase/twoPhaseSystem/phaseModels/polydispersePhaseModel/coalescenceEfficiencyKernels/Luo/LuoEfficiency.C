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

#include "LuoEfficiency.H"
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
    defineTypeNameAndDebug(Luo, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyKernel,
        Luo,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Luo::Luo
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& continuousPhase
)
:
    coalescenceEfficiencyKernel(dict, mesh, continuousPhase),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    epsilonf_
    (
        IOobject
        (
            "Luo:epsilonf",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", sqr(dimVelocity)/dimTime, 0.0)
    ),
    Cvm_
    (
        IOobject
        (
            "Luo:Cvm",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Luo::~Luo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Luo::update
(
    const fluidThermo& thermo,
    const turbulenceModel& turb
)
{
    epsilonf_ = turb.epsilon();
    epsilonf_.max(SMALL);

    const virtualMassModel& virtualMass =
        fluid_.virtualMass(fluid_.phase1());
    Cvm_ = virtualMass.Cvm();
}


Foam::scalar Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Luo::Pc
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli
) const
{
    scalar rhob = fluid_.phase1().rho()[celli];
    scalar rhof = fluid_.phase2().rho()[celli];
    scalar sigma = fluid_.sigma().value();

    scalar xi = min(d1, d2)/max(d1, d2);
    scalar uRel =
        2.0*cbrt(epsilonf_[celli])
       *sqrt(pow(d1, 2.0/3.0) + pow(d2, 2.0/3.0));
    scalar We = min(d1, d2)*rhof*sqr(uRel)/sigma;

    return
        Foam::exp
        (
          - 0.75*sqrt(1.0 + sqr(xi))*sqrt(1.0 + pow3(xi))
           /(sqrt(rhob/rhof + Cvm_[celli])*pow3(1.0 + xi))
           *sqrt(We)
        );
}

// ************************************************************************* //
