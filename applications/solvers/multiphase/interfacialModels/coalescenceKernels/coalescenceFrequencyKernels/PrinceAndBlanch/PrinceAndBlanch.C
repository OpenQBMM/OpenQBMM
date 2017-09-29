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
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceFrequencyKernels
{
    defineTypeNameAndDebug(PrinceAndBlanch, 0);

    addToRunTimeSelectionTable
    (
        coalescenceFrequencyKernel,
        PrinceAndBlanch,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyKernels::PrinceAndBlanch::PrinceAndBlanch
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalescenceFrequencyKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    turbulent_(dict.lookupOrDefault("turbulentCoalescence", false)),
    buoyant_(dict.lookupOrDefault("buoyantCoalescence", true)),
    LS_(dict.lookupOrDefault("laminarShearCoalescence", true))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceFrequencyKernels::PrinceAndBlanch::~PrinceAndBlanch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceFrequencyKernels::PrinceAndBlanch::omega
(
    const label nodei,
    const label nodej
) const
{
    tmp<volScalarField> tmpFreqSrc
    (
        new volScalarField
        (
            IOobject
            (
                "freqSrc",
                fluid_.mesh().time().timeName(),
                fluid_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluid_.mesh(),
            dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0.0)
        )
    );
    volScalarField& freqSrc = tmpFreqSrc.ref();

    const volScalarField& d1 = fluid_.phase1().ds(nodei);
    const volScalarField& d2 = fluid_.phase1().ds(nodej);
    const volScalarField& rho = fluid_.phase2().rho();
    const dimensionedScalar& sigma = fluid_.sigma();
    dimensionedScalar g = mag(fluid_.g());

    if (turbulent_)
    {
        freqSrc == freqSrc
          + 0.089*constant::mathematical::pi*sqr(d1 + d2)
           *sqrt(pow(d1, 2.0/3.0) + pow(d2, 2.0/3.0))
           *cbrt(fluid_.phase2().turbulence().epsilon());
    }
    if (buoyant_)
    {
        freqSrc == freqSrc
          + constant::mathematical::pi/4.0*sqr(d1 + d2)
           *(
               sqrt(2.14*sigma/(d1*rho) + 0.5*g*d1)
             - sqrt(2.14*sigma/(d2*rho) + 0.5*g*d2)
            );
    }
    if (LS_)
    {
        freqSrc == freqSrc
          + 1.0/6.0*pow3(d1 + d2)
           *mag
            (
                fvc::grad(fluid_.phase1().Us(nodei) - fluid_.phase1().Us(nodej))
            );
    }
    return tmpFreqSrc;
}

// ************************************************************************* //
