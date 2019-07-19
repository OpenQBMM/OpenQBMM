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

#include "LuoSvendsenBubble.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include "phaseModel.H"
#include "PhaseCompressibleTurbulenceModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(LuoSvendsenBubble, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        LuoSvendsenBubble,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::LuoSvendsenBubble::
LuoSvendsenBubble
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    breakupKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    rhof_(fluid_.phase2().rho()),
    alphaf_(fluid_.phase2()),
    sigma_(fluid_.sigma()),
    Cf_
    (
        dict.lookupOrDefault
        (
            "Cf",
            dimensionedScalar("Cf", dimless, 0.26)
        )
    ),
    epsilonf_
    (
        IOobject
        (
            "LuoSvendsenBubble:epsilonf",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", sqr(dimVelocity)/dimTime, 0.0)
    ),
    de_
    (
        IOobject
        (
            "LuoSvendsenBubble:de",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::LuoSvendsenBubble::
~LuoSvendsenBubble()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::breakupKernels::LuoSvendsenBubble::
preUpdate()
{
    const phaseModel& phase(fluid_.phase2());
    epsilonf_ = phase.turbulence().epsilon();
    epsilonf_.max(SMALL);

    de_ = pow(pow3(phase.turbulence().nut())/epsilonf_, 0.25);
}


Foam::scalar
Foam::populationBalanceSubModels::breakupKernels::LuoSvendsenBubble::Kb
(
    const scalar& d,
    const label celli,
    const label environment
) const
{
    scalar xi = max(de_[celli]/d, 20.0);

    return
        0.923*alphaf_[celli]*cbrt(epsilonf_[celli]*d)
       *sqr(1.0 + xi)/(sqr(d)*pow(xi, 11.0/3.0))
       *exp
        (
            -12.0*Cf_.value()*sigma_.value()
           /(
               2.045*rhof_[celli]
              *pow(xi, 11.0/3.0)*pow(d, 5.0/3.0)*pow(epsilonf_[celli], 2.0/3.0)
            )
        );
}

// ************************************************************************* //
