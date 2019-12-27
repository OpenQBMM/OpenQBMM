/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 Alberto Passalacqua
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

#include "LuoFrequency.H"
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
namespace coalescenceFrequencyKernels
{
    defineTypeNameAndDebug(Luo, 0);

    addToRunTimeSelectionTable
    (
        coalescenceFrequencyKernel,
        Luo,
        dictionary
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::
coalescenceFrequencyKernels::Luo::Luo
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& continuousPhase
)
:
    coalescenceFrequencyKernel(dict, mesh, continuousPhase),
    epsilonf_
    (
        IOobject
        (
            "Luo:epsilonf",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", sqr(dimVelocity)/dimTime, Zero)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::
coalescenceFrequencyKernels::Luo::~Luo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::populationBalanceSubModels::aggregationKernels::
coalescenceFrequencyKernels::Luo::update
(
    const fluidThermo& thermo,
    const turbulenceModel& turb
)
{
    epsilonf_ = turb.epsilon();
    epsilonf_.max(SMALL);
}


Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::
coalescenceFrequencyKernels::Luo::omega
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli
) const
{
    scalar uRel =
        2.0*cbrt(epsilonf_[celli])
       *sqrt(pow(d1, 2.0/3.0) + pow(d2, 2.0/3.0));

    return constant::mathematical::pi/4.0*sqr(d1 + d2)*uRel;
}

// ************************************************************************* //
