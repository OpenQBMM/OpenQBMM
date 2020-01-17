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

#include "basicVelocityFieldMomentInversion.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicVelocityFieldMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        fieldMomentInversion,
        basicVelocityFieldMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicVelocityFieldMomentInversion::basicVelocityFieldMomentInversion
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
    momentInverter_
    (
        multivariateMomentInversion::New
        (
            dict.subDict("basicVelocityMomentInversion"),
            momentOrders,
            nodeIndexes,
            velocityIndexes
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicVelocityFieldMomentInversion::~basicVelocityFieldMomentInversion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::basicVelocityFieldMomentInversion::invert
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    NotImplemented;
}

void Foam::basicVelocityFieldMomentInversion::invertBoundaryMoments
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    NotImplemented;
}

bool Foam::basicVelocityFieldMomentInversion::invertLocalMoments
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes,
    const label celli,
    const bool fatalErrorOnFailedRealizabilityTest
)
{
    NotImplemented;

    return true;
}

void Foam::basicVelocityFieldMomentInversion::invert
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes
)
{
    const volScalarField& m0(moments(0));

    forAll(m0, celli)
    {
        invertLocalMoments(moments, nodes, celli);
    }

    invertBoundaryMoments(moments, nodes);
}

void Foam::basicVelocityFieldMomentInversion::invertBoundaryMoments
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes
)
{
    // Recover reference to boundaryField of zero-order moment.
    const volScalarField::Boundary& bf = moments[0].boundaryField();

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            multivariateMomentSet momentsToInvert
            (
                moments.size(),
                momentOrders_,
                moments.support()
            );

            // Copying moments from a face
            forAll(momentsToInvert, momenti)
            {
                const labelList& momentOrder = momentOrders_[momenti];

                momentsToInvert(momentOrder)
                        = moments(momentOrder).boundaryField()[patchi][facei];
            }

            // Find quadrature
            momentInverter_().invert(momentsToInvert);

            const mappedList<scalar>& weights(momentInverter_->weights());

            const mappedList<scalarList>& abscissae
            (
                momentInverter_->abscissae()
            );

            const mappedList<vector>& velocityAbscissae
            (
                momentInverter_->velocityAbscissae()
            );

            // Copy quadrature data to boundary face
            forAll(weights, nodei)
            {
                const labelList& nodeIndex = nodeIndexes_[nodei];
                volVelocityNode& node = nodes[nodei];

                volScalarField::Boundary& weightBf =
                    node.primaryWeight().boundaryFieldRef();

                volVectorField::Boundary& velocityAbscissaBf =
                    node.velocityAbscissae().boundaryFieldRef();

                weightBf[patchi][facei] = weights(nodeIndex);

                velocityAbscissaBf[patchi][facei] =
                    velocityAbscissae(nodeIndex);

                forAll(node.scalarIndexes(), cmpt)
                {
                    volScalarField::Boundary& abscissaBf =
                        node.primaryAbscissae()[cmpt].boundaryFieldRef();

                    abscissaBf[patchi][facei] = abscissae(nodeIndex)[cmpt];
                }
            }
        }
    }
}

bool Foam::basicVelocityFieldMomentInversion::invertLocalMoments
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes,
    const label celli,
    const bool fatalErrorOnFailedRealizabilityTest
)
{
    multivariateMomentSet momentsToInvert
    (
        moments.size(),
        momentOrders_,
        moments.support()
    );

    // Copying moments from cell
    forAll(momentsToInvert, momenti)
    {
        const labelList& momentOrder = momentOrders_[momenti];
        momentsToInvert(momentOrder) = moments(momentOrder)[celli];
    }

    if (!momentInverter_().invert(momentsToInvert))
    {
        return false;
    }

    // Recovering quadrature
    const mappedScalarList& weights(momentInverter_().weights());
    const mappedList<scalarList>& abscissae(momentInverter_().abscissae());
    const mappedVectorList& velocityAbscissae
    (
        momentInverter_().velocityAbscissae()
    );

    forAll(nodes, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        volVelocityNode& node(nodes[nodei]);

        node.primaryWeight()[celli] = weights(nodeIndex);
        node.velocityAbscissae()[celli] = velocityAbscissae(nodeIndex);

        forAll(node.scalarIndexes(), cmpt)
        {
            node.primaryAbscissae()[cmpt][celli] = abscissae(nodeIndex)[cmpt];
        }
    }

    return true;
}

// ************************************************************************* //
