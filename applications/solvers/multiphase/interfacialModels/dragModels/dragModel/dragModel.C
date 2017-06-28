/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
2017-05-24 Jeff Heylmun:    Added return functions for acceleration
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

#include "dragModel.H"
#include "phasePair.H"
#include "swarmCorrection.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dragModel, 0);
    defineRunTimeSelectionTable(dragModel, dictionary);
}

const Foam::dimensionSet Foam::dragModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModel::dragModel
(
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair)
{}


Foam::dragModel::dragModel
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair),
    swarmCorrection_
    (
        swarmCorrection::New
        (
            dict.subDict("swarmCorrection"),
            pair
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModel::~dragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModel::CdRe() const
{
    label nNodesi = pair_.dispersed().nNodes();
    label nNodesj = pair_.continuous().nNodes();

    tmp<volScalarField> tCdRe;
    tCdRe =
        CdRe(0, 0)
       *max
        (
            pair_.dispersed().alphas(0)
           *pair_.continuous().alphas(0),
            pair_.dispersed().residualAlpha()
        );

    for (label nodei = 1; nodei < nNodesi; nodei++)
    {
        for (label nodej = 1; nodej < nNodesj; nodej++)
        {
            tCdRe.ref() +=
                CdRe(nodei, nodej)
               *max
                (
                    pair_.dispersed().alphas(0)
                   *pair_.continuous().alphas(0),
                    pair_.dispersed().residualAlpha()
                );
        }
    }
    return
        tCdRe
       /max
        (
            pair_.dispersed()*pair_.continuous(),
            pair_.dispersed().residualAlpha()
        );
}


Foam::tmp<Foam::volScalarField> Foam::dragModel::Ki
(
    const label nodei,
    const label nodej
) const
{
    return
        0.75
       *CdRe(nodei,nodej)
       *swarmCorrection_->Cs(nodei,nodej)
       *pair_.continuous().rho()
       *pair_.continuous().nu()
       /sqr(pair_.dispersed().ds(nodei));
}


Foam::tmp<Foam::volScalarField> Foam::dragModel::Ki() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());
    label nNodesi = pair_.dispersed().nNodes();
    label nNodesj = pair_.continuous().nNodes();

    tmp<volScalarField> tKdi
    (
        new volScalarField
        (
            IOobject
            (
                "totalKdi",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("0", dimK, 0.0)
        )
    );

    for (label nodei = 0; nodei < nNodesi; nodei++)
    {
        for (label nodej = 0; nodej < nNodesj; nodej++)
        {
            tKdi.ref() +=
                K(nodei, nodej)
               *max(pair_.continuous(), pair_.continuous().residualAlpha());
        }
    }
    tKdi.ref() /= max
    (
        pair_.dispersed()*pair_.continuous(),
        pair_.dispersed().residualAlpha()
    );
    return tKdi;
}


Foam::tmp<Foam::volScalarField> Foam::dragModel::K
(
    const label nodei,
    const label nodej
) const
{
    return
        max
        (
            pair_.dispersed().alphas(nodei),
            pair_.dispersed().residualAlpha()
        )*Ki(nodei,nodej);
}


Foam::tmp<Foam::volScalarField> Foam::dragModel::K() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());
    label nNodesi = pair_.dispersed().nNodes();
    label nNodesj = pair_.continuous().nNodes();

    tmp<volScalarField> tKd
    (
        new volScalarField
        (
            IOobject
            (
                "totalKd",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("0", dimK, 0.0)
        )
    );

    for (label nodei = 0; nodei < nNodesi; nodei++)
    {
        for (label nodej = 0; nodej < nNodesj; nodej++)
        {
            tKd.ref() += K(nodei, nodej);
        }
    }
    return tKd;
}


Foam::tmp<Foam::surfaceScalarField> Foam::dragModel::Kf
(
    const label nodei,
    const label nodej
) const
{
    return
        max
        (
            fvc::interpolate(pair_.dispersed().alphas(nodei)),
            pair_.dispersed().residualAlpha()
        )*fvc::interpolate(Ki(nodei,nodej));
}


bool Foam::dragModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
