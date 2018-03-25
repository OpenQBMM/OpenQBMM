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

#include "hyperbolicFieldMomentInversion.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hyperbolicFieldMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        fieldMomentInversion,
        hyperbolicFieldMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hyperbolicFieldMomentInversion::hyperbolicFieldMomentInversion
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
        new hyperbolicConditionalMomentInversion
        (
            dict.subDict("basicMomentInversion"),
            momentOrders[0].size()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hyperbolicFieldMomentInversion::~hyperbolicFieldMomentInversion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hyperbolicFieldMomentInversion::invert
(
    const volUnivariateMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    NotImplemented;
}

void Foam::hyperbolicFieldMomentInversion::invertBoundaryMoments
(
    const volUnivariateMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    NotImplemented;
}

bool Foam::hyperbolicFieldMomentInversion::invertLocalMoments
(
    const volUnivariateMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes,
    const label celli,
    const bool fatalErrorOnFailedRealizabilityTest
)
{
    NotImplemented;

    return true;
}

void Foam::hyperbolicFieldMomentInversion::invert
(
    const volVectorMomentFieldSet& moments,
    mappedPtrList<volVectorNode>& nodes
)
{
    const volVectorMoment& m0(moments[0]);

    forAll(m0, celli)
    {
        invertLocalMoments(moments, nodes, celli);
    }

    invertBoundaryMoments(moments, nodes);
}

void Foam::hyperbolicFieldMomentInversion::invertBoundaryMoments
(
    const volVectorMomentFieldSet& moments,
    mappedPtrList<volVectorNode>& nodes
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
            const mappedList<vector>& abscissae(momentInverter_->abscissae());

            // Copy quadrature data to boundary face
            forAll(weights, nodei)
            {
                const labelList& nodeIndex = nodeIndexes_[nodei];
                volVectorNode& node = nodes[nodei];

                volScalarField::Boundary& weightBf
                        = node.primaryWeight().boundaryFieldRef();

                volVectorField::Boundary& abscissaBf
                        = node.primaryAbscissa().boundaryFieldRef();

                weightBf[patchi][facei] = weights(nodeIndex);
                abscissaBf[patchi][facei] = abscissae(nodeIndex);
            }
        }
    }
}

bool Foam::hyperbolicFieldMomentInversion::invertLocalMoments
(
    const volVectorMomentFieldSet& moments,
    mappedPtrList<volVectorNode>& nodes,
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
    const mappedList<scalar>& weights(momentInverter_().weights());
    const mappedList<vector>& abscissae(momentInverter_().abscissae());

    forAll(weights, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        volVectorNode& node(nodes[nodei]);

        node.primaryWeight()[celli] = weights(nodeIndex);
        node.primaryAbscissa()[celli] = abscissae(nodeIndex);
    }

    return true;
}

// ************************************************************************* //
