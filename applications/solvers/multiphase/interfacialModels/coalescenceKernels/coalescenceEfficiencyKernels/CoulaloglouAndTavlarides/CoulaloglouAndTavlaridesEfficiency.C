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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::
CoulaloglouAndTavlarides
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalescenceEfficiencyKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    Ceff_(dict.lookup("Ceff")),
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

Foam::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::
~CoulaloglouAndTavlarides()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::update()
{
    const phaseModel& phase(fluid_.phase1());
    epsilonf_ = phase.turbulence().epsilon();
    epsilonf_.max(SMALL);
    muf_ = fluid_.phase2().mu();
}


Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::Pc
(
    const label nodei,
    const label nodej
) const
{
    scalar pi = constant::mathematical::pi;
    volScalarField v1(pow3(fluid_.phase1().ds(nodei))*pi/6.0);
    volScalarField v2(pow3(fluid_.phase1().ds(nodej))*pi/6.0);
    const volScalarField& rho = fluid_.phase2().rho();
    const dimensionedScalar& sigma = fluid_.sigma();

    return
        Foam::exp
        (
          - Ceff_
           *muf_*epsilonf_*rho/sqr(sigma)
           *pow4(cbrt(v1)*cbrt(v2)/(cbrt(v1) + cbrt(v2)))
        );
}


Foam::scalar Foam::coalescenceEfficiencyKernels::CoulaloglouAndTavlarides::Pc
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
