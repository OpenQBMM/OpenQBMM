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

#include "LuoSvendsen.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include "phaseModel.H"
#include "PhaseCompressibleTurbulenceModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace bubbleBreakupKernels
{
    defineTypeNameAndDebug(LuoSvendsen, 0);

    addToRunTimeSelectionTable
    (
        bubbleBreakupKernel,
        LuoSvendsen,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubbleBreakupKernels::LuoSvendsen::LuoSvendsen
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    bubbleBreakupKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    C1_
    (
        dict.lookupOrDefault
        (
            "C1",
            dimensionedScalar("C1", dimless, 0.04)
        )
    ),
    C2_
    (
        dict.lookupOrDefault
        (
            "C2",
            dimensionedScalar("C2", dimless, 0.01)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bubbleBreakupKernels::LuoSvendsen::~LuoSvendsen()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::bubbleBreakupKernels::LuoSvendsen::Kb(const label nodei) const
{
    volScalarField epsilon
    (
        "epsilon",
        fluid_.phase2().turbulence().epsilon()
    );
    epsilon.max(SMALL);

    const volScalarField& d = fluid_.phase1().ds(nodei);
    const volScalarField& rho1 = fluid_.phase1().rho();
    const volScalarField& rho2 = fluid_.phase2().rho();
    const volScalarField& mu = fluid_.phase2().mu();
    const dimensionedScalar& sigma = fluid_.sigma();

    tmp<volScalarField> breakupSource =
        Cb_*cbrt(epsilon)
       *Foam::erfc
        (
            Foam::sqrt
            (
                C1_*sigma
               /(rho2*pow(epsilon, 2.0/3.0)*pow(d, 5.0/3.0))
              + C2_*mu
               /(sqrt(rho1*rho2)*cbrt(epsilon)*pow(d, 4.0/3.0))
            )
        );
    breakupSource.ref().dimensions().reset(inv(dimTime));
    return breakupSource;
}

// ************************************************************************* //
