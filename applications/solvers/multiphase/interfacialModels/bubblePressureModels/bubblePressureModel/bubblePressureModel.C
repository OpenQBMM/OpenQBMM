/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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

#include "bubblePressureModel.H"
#include "phasePair.H"
#include "fvc.H"
#include "fvcFlux.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "virtualMassModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bubblePressureModel, 0);
    defineRunTimeSelectionTable(bubblePressureModel, dictionary);
}

const Foam::dimensionSet Foam::bubblePressureModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bubblePressureModel::bubblePressureModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair),
    Cbp_
    (
        dimensionedScalar::lookupOrDefault
        (
            "Cbp",
            dict,
            dimless,
            1.0
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::bubblePressureModel::~bubblePressureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::bubblePressureModel::nuEff
(
    const label nodei,
    const label nodej
) const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    const virtualMassModel&
        virtualMass
        (
            mesh.lookupObject<virtualMassModel>
            (
                IOobject::groupName(virtualMassModel::typeName, pair_.name())
            )
        );

    return
        pair_.dispersed().alphas(nodei)
       /max
        (
            pair_.continuous().alphas(nodej),
            pair_.continuous().residualAlpha()
        )
       *(
            pair_.dispersed().rho()/pair_.continuous().rho()
          + virtualMass.Cvm(nodei, nodej)
        )*pair_.dispersed().nu();
}

Foam::tmp<Foam::volVectorField> Foam::bubblePressureModel::divDevRhoReff
(
    const label nodei,
    const label nodej
) const
{
    return
        fvc::div
        (
            pair_.continuous().rho()
           *pair_.continuous().alphas(nodej)
           *nuEff(nodei, nodej)
           *dev(twoSymm(fvc::grad(pair_.continuous().Us(nodej))))
        );
}


Foam::tmp<Foam::volVectorField> Foam::bubblePressureModel::Fi
(
    const label nodei,
    const label nodej
) const
{
    return Cbp_*(divDevRhoReff(nodei, nodej) - fvc::grad(pb(nodei, nodej)));
}


Foam::tmp<Foam::volVectorField> Foam::bubblePressureModel::F
(
    const label nodei,
    const label nodej
) const
{
    return pair_.dispersed().alphas(nodei)*Fi(nodei, nodej);
}


Foam::tmp<Foam::surfaceScalarField> Foam::bubblePressureModel::Ff
(
    const label nodei,
    const label nodej
) const
{
    return
        fvc::interpolate(pair_.dispersed().alphas(nodei))
       *fvc::flux(Fi(nodei, nodej));
}


// ************************************************************************* //
