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

#include "LuoFrequency.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyKernels::Luo::Luo
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
            "Luo:epsilonf",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", sqr(dimVelocity)/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyKernels::Luo::~Luo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coalescenceFrequencyKernels::Luo::update()
{
    const phaseModel& phase(fluid_.phase1());
    epsilonf_ = phase.turbulence().epsilon();
    epsilonf_.max(SMALL);
}

Foam::tmp<Foam::volScalarField>
Foam::coalescenceFrequencyKernels::Luo::omega
(
    const label nodei,
    const label nodej
) const
{
    const volScalarField& d1 = fluid_.phase1().ds(nodei);
    const volScalarField& d2 = fluid_.phase1().ds(nodej);
    const volScalarField& rho = fluid_.phase2().rho();
    const dimensionedScalar& sigma = fluid_.sigma();

    volScalarField xi("xi", min(d1, d2)/max(d1, d2));
    volScalarField uRel
    (
        "uRel",
        2.0*cbrt(epsilonf_)*sqrt(pow(d1, 2.0/3.0) + pow(d2, 2.0/3.0))
    );

    return constant::mathematical::pi/4.0*sqr(d1 + d2)*uRel;
}

Foam::scalar
Foam::coalescenceFrequencyKernels::Luo::omega
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    scalar d1 = fluid_.phase1().ds(nodei)[celli];
    scalar d2 = fluid_.phase1().ds(nodej)[celli];

    scalar xi = min(d1, d2)/max(d1, d2);
    scalar uRel =
        2.0*cbrt(epsilonf_[celli])*sqrt(pow(d1, 2.0/3.0) + pow(d2, 2.0/3.0));

    return constant::mathematical::pi/4.0*sqr(d1 + d2)*uRel;
}

// ************************************************************************* //
