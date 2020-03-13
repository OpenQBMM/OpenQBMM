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

#include "ChestersEfficiency.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace coalescenceEfficiencyKernels
{
    defineTypeNameAndDebug(Chesters, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyKernel,
        Chesters,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Chesters::Chesters
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& continuousPhase
)
:
    coalescenceEfficiencyKernel(dict, mesh, continuousPhase),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    Ceff_("Ceff", dimless, dict),
    ReExp_("ReExp", dimless, dict),
    WeExp_("WeExp", dimless, dict),
    theta_
    (
        IOobject
        (
            "Chesters:theta",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Chesters::~Chesters()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Chesters::update
(
    const fluidThermo& thermo,
    const turbulenceModel& turb
)
{
    const phasePair& pair = fluid_.pair1In2();
    theta_ =
        Ceff_
       *pow(max(pair.Re(), SMALL), ReExp_)
       *pow(max(pair.We(), SMALL), WeExp_);
}


Foam::scalar Foam::populationBalanceSubModels::aggregationKernels::coalescenceEfficiencyKernels::Chesters::Pc
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli
) const
{
    const phaseModel& phase1 = fluid_.phase1();
    const phaseModel& phase2 = fluid_.phase2();

    scalar Weij
    (
        phase2.rho()[celli]*d1*magSqr(Ur)/fluid_.sigma().value()
    );
    scalar xi(d1/d2);

    return
        Foam::exp
        (
          - theta_[celli]*sqrt(Weij)
           *sqrt(0.75*(1.0 + sqr(xi))*(1.0 + pow3(xi)))
           /(
                sqrt
                (
                    phase1.rho()[celli]/phase2.rho()[celli]
                  + 0.5
                )
               *pow3(1.0 + xi)
            )
        );
}

// ************************************************************************* //
