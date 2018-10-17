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

#include "CoulaloglouAndTavlaridesFrequency.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
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
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::
CoulaloglouAndTavlarides
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalescenceFrequencyKernel(dict, mesh),
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

Foam::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::
~CoulaloglouAndTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::update()
{
    const phaseModel& phase(fluid_.phase1());
    epsilonf_ = phase.turbulence().epsilon();
    epsilonf_.max(SMALL);
}


Foam::tmp<Foam::volScalarField>
Foam::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::omega
(
    const label nodei,
    const label nodej
) const
{
    scalar pi = constant::mathematical::pi;
    volScalarField v1(pow3(fluid_.phase1().ds(nodei))*pi/6.0);
    volScalarField v2(pow3(fluid_.phase1().ds(nodej))*pi/6.0);

    return
        (pow(v1, 2.0/3.0) + pow(v2, 2.0/3.0))
       *sqrt(pow(v1, 2.0/9.0) + pow(v2, 2.0/9.0))
       *cbrt(epsilonf_)/(1.0 + fluid_.phase1());
}


Foam::scalar Foam::coalescenceFrequencyKernels::CoulaloglouAndTavlarides::omega
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    scalar d1 = fluid_.phase1().ds(nodei)[celli];
    scalar d2 = fluid_.phase1().ds(nodej)[celli];

    scalar pi = constant::mathematical::pi;
    scalar v1 = pow3(d1)*pi/6.0;
    scalar v2 = pow3(d2)*pi/6.0;

    return
        (pow(v1, 2.0/3.0) + pow(v2, 2.0/3.0))
       *sqrt(pow(v1, 2.0/9.0) + pow(v2, 2.0/9.0))
       *cbrt(epsilonf_[celli])/(1.0 + fluid_.phase1()[celli]);
}

// ************************************************************************* //
