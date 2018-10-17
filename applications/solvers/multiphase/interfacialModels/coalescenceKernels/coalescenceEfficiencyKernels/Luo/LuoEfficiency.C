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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyKernels::Luo::Luo
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalescenceEfficiencyKernel(dict, mesh),
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

Foam::coalescenceEfficiencyKernels::Luo::~Luo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coalescenceEfficiencyKernels::Luo::update()
{
    const phaseModel& phase(fluid_.phase1());
    epsilonf_ = phase.turbulence().epsilon();
    epsilonf_.max(SMALL);

    const virtualMassModel& virtualMass = fluid_.virtualMass(fluid_.phase1());
    Cvm_ = virtualMass.Cvm();
}


Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyKernels::Luo::Pc
(
    const label nodei,
    const label nodej
) const
{
    const volScalarField& d1 = fluid_.phase1().ds(nodei);
    const volScalarField& d2 = fluid_.phase1().ds(nodej);
    const volScalarField& rhob = fluid_.phase1().rho();
    const volScalarField& rhof = fluid_.phase2().rho();
    const dimensionedScalar& sigma = fluid_.sigma();

    volScalarField xi("xi", min(d1, d2)/max(d1, d2));
    volScalarField uRel
    (
        "uRel",
        2.0*cbrt(epsilonf_)*sqrt(pow(d1, 2.0/3.0) + pow(d2, 2.0/3.0))
    );
    volScalarField We
    (
        "We",
        min(d1, d2)*rhof*sqr(uRel)/sigma
    );

    return
        Foam::exp
        (
          - 0.75*sqrt(1.0 + sqr(xi))*sqrt(1.0 + pow3(xi))
           /(sqrt(rhob/rhof + Cvm_)*pow3(1.0 + xi))
           *sqrt(We)
        );
}


Foam::scalar Foam::coalescenceEfficiencyKernels::Luo::Pc
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    scalar d1 = fluid_.phase1().ds(nodei)[celli];
    scalar d2 = fluid_.phase1().ds(nodej)[celli];
    scalar rhob = fluid_.phase1().rho()[celli];
    scalar rhof = fluid_.phase2().rho()[celli];
    scalar sigma = fluid_.sigma().value();

    scalar xi = min(d1, d2)/max(d1, d2);
    scalar uRel =
        2.0*cbrt(epsilonf_[celli])*sqrt(pow(d1, 2.0/3.0) + pow(d2, 2.0/3.0));
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
