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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyKernels::Chesters::Chesters
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalescenceEfficiencyKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    Ceff_(dict.lookup("Ceff")),
    ReExp_(dict.lookup("ReExp")),
    WeExp_(dict.lookup("WeExp")),
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

Foam::coalescenceEfficiencyKernels::Chesters::~Chesters()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coalescenceEfficiencyKernels::Chesters::update()
{
    const phasePair& pair = fluid_.pair1In2();
    theta_ =
        Ceff_
       *pow(max(pair.Re(), SMALL), ReExp_)
       *pow(max(pair.We(), SMALL), WeExp_);
}

Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyKernels::Chesters::Pc
(
    const label nodei,
    const label nodej
) const
{

    const phaseModel& phase1 = fluid_.phase1();
    const phaseModel& phase2 = fluid_.phase2();

    const volScalarField& di = fluid_.phase1().ds(nodei);
    const volScalarField& dj = fluid_.phase1().ds(nodej);

    volScalarField Weij
    (
        phase2.rho()
       *di
       *magSqr(phase1.Us(nodei) - phase1.Us(nodej))
       /fluid_.sigma()
    );
    volScalarField xi(di/dj);

    return
        Foam::exp
        (
          - theta_*sqrt(Weij)
           *sqrt(0.75*(1.0 + sqr(xi))*(1.0 + pow3(xi)))
           /(
                sqrt(fluid_.phase1().rho()/fluid_.phase2().rho() + 0.5)
               *pow3(1.0 + xi)
            )
        );
}


Foam::scalar Foam::coalescenceEfficiencyKernels::Chesters::Pc
(
    const label celli,
    const label nodei,
    const label nodej
) const
{
    const phaseModel& phase1 = fluid_.phase1();
    const phaseModel& phase2 = fluid_.phase2();

    scalar di = fluid_.phase1().ds(nodei)[celli];
    scalar dj = fluid_.phase1().ds(nodej)[celli];

    scalar Weij
    (
        phase2.rho()[celli]
       *di
       *magSqr(phase1.Us(nodei)[celli] - phase1.Us(nodej)[celli])
       /fluid_.sigma().value()
    );
    scalar xi(di/dj);

    return
        Foam::exp
        (
          - theta_[celli]*sqrt(Weij)
           *sqrt(0.75*(1.0 + sqr(xi))*(1.0 + pow3(xi)))
           /(
                sqrt
                (
                    fluid_.phase1().rho()[celli]/fluid_.phase2().rho()[celli]
                  + 0.5
                )
               *pow3(1.0 + xi)
            )
        );
}

// ************************************************************************* //
