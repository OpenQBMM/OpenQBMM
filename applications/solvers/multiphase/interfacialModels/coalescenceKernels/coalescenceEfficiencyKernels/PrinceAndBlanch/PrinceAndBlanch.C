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

#include "PrinceAndBlanch.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceEfficiencyKernels
{
    defineTypeNameAndDebug(PrinceAndBlanch, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyKernel,
        PrinceAndBlanch,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyKernels::PrinceAndBlanch::
PrinceAndBlanch
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalescenceEfficiencyKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    ho_
    (
        dimensionedScalar::lookupOrDefault
        (
            "h0",
            dict,
            dimLength,
            1e-3
        )
    ),
    hf_
    (
        dimensionedScalar::lookupOrDefault
        (
            "hf",
            dict,
            dimLength,
            1e-6
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyKernels::PrinceAndBlanch::
~PrinceAndBlanch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyKernels::PrinceAndBlanch::Pc
(
    const label nodei,
    const label nodej
) const
{
    const volScalarField& d1 = fluid_.phase1().ds(nodei);
    const volScalarField& d2 = fluid_.phase1().ds(nodej);
    const volScalarField& rho = fluid_.phase2().rho();
    const volScalarField& epsilon = fluid_.phase2().turbulence().epsilon();
    const dimensionedScalar& sigma = fluid_.sigma();

    volScalarField rij("rij", 0.5/(2.0/d1 + 2.0/d2));

    return
        Foam::exp
        (
          - sqrt(rho*pow3(rij)/(16.0*sigma))*log(ho_/hf_)
           /(pow(rij, 2.0/3.0)/pow(epsilon, 1.0/3.0))
        );
}

// ************************************************************************* //
