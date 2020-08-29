/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2020 Alberto Passalacqua - Implemented Cunningham correction.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "CunninghamSchillerNaumann.H"
#include "physicoChemicalConstants.H"
#include "phasePair.H"
#include "SchillerNaumann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(CunninghamSchillerNaumann, 0);
    addToRunTimeSelectionTable(dragModel, CunninghamSchillerNaumann, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::CunninghamSchillerNaumann::CunninghamSchillerNaumann
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    SchillerNaumann_
    (
        new SchillerNaumann
        (
            dict,
            pair,
            false
        )
    ),
    residualRe_("residualRe", dimless, dict),
    A1_(dict.lookupOrDefault<scalar>("A1", 1.257)),
    A2_(dict.lookupOrDefault<scalar>("A2", 0.4)),
    A3_(dict.lookupOrDefault<scalar>("A3", 0.55)),
    M_(dict.lookupOrDefault<scalar>("M", 15.9994))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::CunninghamSchillerNaumann::~CunninghamSchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::CunninghamSchillerNaumann::CdRe
(
    const label nodei,
    const label nodej
) const
{
    const volScalarField& dp = pair_.dispersed().d(nodei);
    const volScalarField& mu = pair_.continuous().thermo().mu();
    const volScalarField& T = pair_.continuous().thermo().T();
    const volScalarField& p = pair_.continuous().thermo().p();

    // Calculate the mean free path using the kinetic theory of gases made
    // of hard spheres. The gas is assumed to have the same viscosity and the
    // same molecular weight of the real gas.
    volScalarField meanFreePath
    (
        mu*sqrt
        (
            Foam::constant::mathematical::pi
            *Foam::constant::physicoChemical::R*T/(2.0*M_)
        )/p
    );

    // Calculate Cunnigham correction factor to account for rarefaction effects
    volScalarField CunninghamCorrection
    (
        scalar(1) + 2.0*meanFreePath*(A1_ + A2_*exp(-A3_*dp/meanFreePath))/dp
    );

    // Return the standard Schiller-Naumann coefficient divided by the
    // Cunnigham correction factor.
    return SchillerNaumann_->CdRe(nodei, nodej)/CunninghamCorrection;
}


// ************************************************************************* //
