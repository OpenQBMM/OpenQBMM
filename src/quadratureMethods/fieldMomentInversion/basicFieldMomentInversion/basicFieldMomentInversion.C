/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
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

#include "basicFieldMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicFieldMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        fieldMomentInversion,
        basicFieldMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicFieldMomentInversion::basicFieldMomentInversion
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes,
    const label nSecondaryNodes
)
:
    fieldMomentInversion
    (
        dict,
        mesh,
        momentOrders,
        nodeIndexes,
        velocityIndexes,
        nSecondaryNodes
    ),
    minKnownAbscissa_(dict.lookupOrDefault<scalar>("minKnownAbscissa", 0)),
    maxKnownAbscissa_(dict.lookupOrDefault<scalar>("maxKnownAbscissa", 1)),
    nFixedQuadraturePoints_(0),
    momentInverter_
    (
        univariateMomentInversion::New(dict.subDict("basicMomentInversion"))
    )
{
    static word inversionType = momentInverter_().type();

    if (inversionType == "GaussRadau")
    {
        nFixedQuadraturePoints_ = 1;
    }
    else if (inversionType == "GaussLobatto")
    {
        nFixedQuadraturePoints_ = 2;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicFieldMomentInversion::~basicFieldMomentInversion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::basicFieldMomentInversion::invert
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    const volScalarField& m0(moments(0));

    forAll(m0, celli)
    {
        invertLocalMoments(moments, nodes, celli);
    }

    invertBoundaryMoments(moments, nodes);
}

void Foam::basicFieldMomentInversion::invertBoundaryMoments
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    // Recover reference to boundaryField of zero-order moment.
    const volScalarField::Boundary& bf = moments[0].boundaryField();

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            univariateMomentSet momentsToInvert
            (
                moments.size(),
                moments.support(),
                scalar(0),                  // Initial value
                nFixedQuadraturePoints_
            );

            // Copying moments from a face
            forAll(momentsToInvert, momenti)
            {
                momentsToInvert[momenti]
                        = moments[momenti].boundaryField()[patchi][facei];
            }

            // Find quadrature
            momentInverter_().invert
            (
                momentsToInvert,
                minKnownAbscissa_,
                maxKnownAbscissa_
            );

            label maxNodes = nodes.size();
            label actualNodes = momentInverter_().nNodes();

            // Copy quadrature data to boundary face
            for (label nodei = 0; nodei < maxNodes; nodei++)
            {
                volScalarNode& node = nodes[nodei];

                volScalarField::Boundary& weightBf
                        = node.primaryWeight().boundaryFieldRef();

                volScalarField::Boundary& abscissaBf
                        = node.primaryAbscissae()[0].boundaryFieldRef();

                if (nodei < actualNodes)
                {
                    weightBf[patchi][facei]
                            = momentInverter_().weights()[nodei];

                    abscissaBf[patchi][facei]
                            = momentInverter_().abscissae()[nodei];
                }
                else
                {
                    weightBf[patchi][facei] = 0.0;
                    abscissaBf[patchi][facei] = 0.0;
                }
            }
        }
    }
}

bool Foam::basicFieldMomentInversion::invertLocalMoments
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes,
    const label celli,
    const bool fatalErrorOnFailedRealizabilityTest
)
{
    univariateMomentSet momentsToInvert
    (
        moments.size(),
        moments.support(),
        scalar(0),                  // Initial value
        nFixedQuadraturePoints_
    );

    // Copying moments from cell
    forAll(momentsToInvert, momenti)
    {
        momentsToInvert[momenti] = moments[momenti][celli];
    }

    if (!fatalErrorOnFailedRealizabilityTest)
    {
        if (!momentsToInvert.isRealizable(fatalErrorOnFailedRealizabilityTest))
        {
            return false;
        }
    }

    // Find quadrature
    momentInverter_().invert
    (
        momentsToInvert,
        minKnownAbscissa_,
        maxKnownAbscissa_
    );

    label maxNodes = nodes.size();
    label actualNodes = momentInverter_().nNodes();

    // Recovering quadrature
    const scalarList& weights(momentInverter_().weights());
    const scalarList& abscissae(momentInverter_().abscissae());

    for (label nodei = 0; nodei < maxNodes; nodei++)
    {
        volScalarNode& node(nodes[nodei]);

        if (nodei < actualNodes)
        {
            node.primaryWeight()[celli] = weights[nodei];
            node.primaryAbscissae()[0][celli] = abscissae[nodei];
        }
        else
        {
            node.primaryWeight()[celli] = 0.0;
            node.primaryAbscissae()[0][celli] = 0.0;
        }
    }

    return true;
}

void Foam::basicFieldMomentInversion::invert
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes
)
{
    NotImplemented;
}

void Foam::basicFieldMomentInversion::invertBoundaryMoments
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes
)
{
    NotImplemented;
}

bool Foam::basicFieldMomentInversion::invertLocalMoments
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes,
    const label celli,
    const bool fatalErrorOnFailedRealizabilityTest
)
{
    NotImplemented;

    return true;
}
// ************************************************************************* //
