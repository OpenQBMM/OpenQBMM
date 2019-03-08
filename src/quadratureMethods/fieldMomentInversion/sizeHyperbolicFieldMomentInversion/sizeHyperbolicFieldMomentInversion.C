/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
     \\/     M anipulation  |
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

#include "sizeHyperbolicFieldMomentInversion.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sizeHyperbolicFieldMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        fieldMomentInversion,
        sizeHyperbolicFieldMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sizeHyperbolicFieldMomentInversion::sizeHyperbolicFieldMomentInversion
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const label nSecondaryNodes
)
:
    fieldMomentInversion
    (
        dict,
        mesh,
        momentOrders,
        nodeIndexes,
        nSecondaryNodes
    ),
    momentInverter_
    (
        new sizeCHyQMOM
        (
            dict.subDict("basicMomentInversion"),
            momentOrders,
            nodeIndexes
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sizeHyperbolicFieldMomentInversion::~sizeHyperbolicFieldMomentInversion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sizeHyperbolicFieldMomentInversion::invert
(
    const volMomentFieldSet& moments,
    mappedPtrList<volNode>& nodes
)
{
    NotImplemented;
}

void Foam::sizeHyperbolicFieldMomentInversion::invertBoundaryMoments
(
    const volMomentFieldSet& moments,
    mappedPtrList<volNode>& nodes
)
{
    NotImplemented;
}

bool Foam::sizeHyperbolicFieldMomentInversion::invertLocalMoments
(
    const volMomentFieldSet& moments,
    mappedPtrList<volNode>& nodes,
    const label celli,
    const bool fatalErrorOnFailedRealizabilityTest
)
{
    NotImplemented;

    return true;
}

void Foam::sizeHyperbolicFieldMomentInversion::invert
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

void Foam::sizeHyperbolicFieldMomentInversion::invertBoundaryMoments
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
            const mappedList<scalar>& sizeAbscissae
            (
                momentInverter_->sizeAbscissae()
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

                volScalarField::Boundary& sizeAbscissaBf =
                    node.primaryAbscissae()[0].boundaryFieldRef();

                volVectorField::Boundary& velocityAbscissaBf =
                    node.velocityAbscissae().boundaryFieldRef();

                weightBf[patchi][facei] = weights(nodeIndex);
                sizeAbscissaBf[patchi][facei] = sizeAbscissae(nodeIndex);
                velocityAbscissaBf[patchi][facei] =
                    velocityAbscissae(nodeIndex);
            }
        }
    }
}

bool Foam::sizeHyperbolicFieldMomentInversion::invertLocalMoments
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

//     if (!fatalErrorOnFailedRealizabilityTest)
//     {
//         if (!momentsToInvert.isRealizable(fatalErrorOnFailedRealizabilityTest))
//         {
//             return false;
//         }
//     }

    // Find quadrature
    momentInverter_().invert(momentsToInvert);

    // Recovering quadrature
    const mappedScalarList& weights(momentInverter_().weights());
    const mappedScalarList& sizeAbscissae(momentInverter_().sizeAbscissae());
    const mappedVectorList& velocityAbscissae
    (
        momentInverter_().velocityAbscissae()
    );

    forAll(weights, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        volVelocityNode& node(nodes[nodei]);

        node.primaryWeight()[celli] = weights(nodeIndex);
        node.primaryAbscissae()[0][celli] = sizeAbscissae(nodeIndex);
        node.velocityAbscissae()[celli] = velocityAbscissae(nodeIndex);
    }

    return true;
}

// ************************************************************************* //
