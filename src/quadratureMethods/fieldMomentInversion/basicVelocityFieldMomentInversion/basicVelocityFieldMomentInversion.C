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
    Copyright (C) 2019-2026 Alberto Passalacqua
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
    const labelList& velocityIndexes
)
:
    fieldMomentInversion
    (
        dict,
        mesh,
        momentOrders,
        nodeIndexes,
        velocityIndexes
    ),
    momentsToInvert_(nullptr),
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

Foam::multivariateMomentSet&
Foam::basicVelocityFieldMomentInversion::localMomentSet
(
    const volVelocityMomentFieldSet& moments
)
{
    if (!momentsToInvert_)
    {
        momentsToInvert_.reset
        (
            new multivariateMomentSet
            (
                moments.size(),
                momentOrders_,
                moments.supports(),
                momentInverter_().smallM0(),
                momentInverter_().smallZeta()
            )
        );
    }
    else if (momentsToInvert_().nMoments() != moments.size())
    {
        FatalErrorInFunction
            << "The moment field set is inconsistent with the moment set "
            << "used to invert it." << nl
            << "    Number of moments: " << moments.size()
            << ", expected " << momentsToInvert_().nMoments() << nl
            << exit(FatalError);
    }

    return momentsToInvert_();
}


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

    label nFailedCells = 0;
    label firstFailedCell = -1;

    forAll(m0, celli)
    {
        if (!invertLocalMoments(moments, nodes, celli))
        {
            if (nFailedCells == 0)
            {
                firstFailedCell = celli;
            }

            nFailedCells++;
        }
    }

    if (nFailedCells > 0)
    {
        WarningInFunction
            << "The moments of " << nFailedCells << " of the " << m0.size()
            << " cells of the mesh could not be inverted, and the "
            << "quadrature of those cells was set to zero." << nl
            << "    First of them: cell " << firstFailedCell
            << ", of zero-order moment " << m0[firstFailedCell] << nl
            << endl;
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

    multivariateMomentSet& momentsToInvert(localMomentSet(moments));

    label nFailedFaces = 0;
    label firstFailedPatch = -1;
    label firstFailedFace = -1;

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            // Copying moments from a face
            forAll(momentsToInvert, momenti)
            {
                const labelList& momentOrder = momentOrders_[momenti];

                momentsToInvert(momentOrder)
                        = moments(momentOrder).boundaryField()[patchi][facei];
            }

            // Find quadrature. A face whose moments have collapsed is
            // emptied and counted, as a cell of the mesh is.
            const bool inverted = momentInverter_().invert(momentsToInvert);

            if (!inverted)
            {
                if (nFailedFaces == 0)
                {
                    firstFailedPatch = patchi;
                    firstFailedFace = facei;
                }

                nFailedFaces++;
            }

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
                    node.weight().boundaryFieldRef();

                volVectorField::Boundary& velocityAbscissaBf =
                    node.velocityAbscissae().boundaryFieldRef();

                weightBf[patchi][facei] =
                    inverted ? weights(nodeIndex) : Zero;

                velocityAbscissaBf[patchi][facei] =
                    inverted ? velocityAbscissae(nodeIndex) : vector::zero;

                forAll(node.scalarIndexes(), cmpt)
                {
                    volScalarField::Boundary& abscissaBf =
                        node.abscissae()[cmpt].boundaryFieldRef();

                    abscissaBf[patchi][facei] =
                        inverted ? abscissae(nodeIndex)[cmpt] : Zero;
                }
            }
        }
    }

    if (nFailedFaces > 0)
    {
        WarningInFunction
            << "The moments of " << nFailedFaces << " faces of the boundary "
            << "could not be inverted, and the quadrature of those faces "
            << "was set to zero." << nl
            << "    First of them: face " << firstFailedFace << " of patch "
            << bf[firstFailedPatch].patch().name() << nl
            << endl;
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
    multivariateMomentSet& momentsToInvert(localMomentSet(moments));

    // Copying moments from cell
    forAll(momentsToInvert, momenti)
    {
        const labelList& momentOrder = momentOrders_[momenti];
        momentsToInvert(momentOrder) = moments(momentOrder)[celli];
    }

    if (!momentInverter_().invert(momentsToInvert))
    {
        // The caller decides what a cell that cannot be inverted means.
        //
        // The adaptive ODE solver asks for the failure to be reported so
        // that it can retry the step, and it retries from what the cell
        // still holds, so its quadrature is left alone.
        if (!fatalErrorOnFailedRealizabilityTest)
        {
            return false;
        }

        // A sweep over the mesh has nothing to retry with. A moment set
        // that no quadrature represents is one that has collapsed, which
        // for the vanishing moment sets this happens on means the cell has
        // been emptied, so the quadrature is set to zero rather than left
        // carrying the one the previous update wrote. The sweep counts the
        // cells and reports them.
        forAll(nodes, nodei)
        {
            volVelocityNode& node(nodes[nodei]);

            node.weight()[celli] = Zero;
            node.velocityAbscissae()[celli] = Zero;

            forAll(node.scalarIndexes(), cmpt)
            {
                node.abscissae()[cmpt][celli] = Zero;
            }
        }

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

        node.weight()[celli] = weights(nodeIndex);
        node.velocityAbscissae()[celli] = velocityAbscissae(nodeIndex);

        forAll(node.scalarIndexes(), cmpt)
        {
            node.abscissae()[cmpt][celli] = abscissae(nodeIndex)[cmpt];
        }
    }

    return true;
}

Foam::scalar Foam::basicVelocityFieldMomentInversion::smallM0() const
{
    return momentInverter_().smallM0();
}

Foam::scalar Foam::basicVelocityFieldMomentInversion::smallZeta() const
{
    return momentInverter_().smallZeta();
}

// ************************************************************************* //
