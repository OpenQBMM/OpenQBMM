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

#include "extendedFieldMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedFieldMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        fieldMomentInversion,
        extendedFieldMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedFieldMomentInversion::extendedFieldMomentInversion
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
        extendedMomentInversion::New
        (
            dict.subDict("extendedMomentInversion"),
            momentOrders.size(),
            nSecondaryNodes
        )
    )
{
    extended_ = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedFieldMomentInversion::~extendedFieldMomentInversion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::extendedFieldMomentInversion::invert
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

void Foam::extendedFieldMomentInversion::invertBoundaryMoments
(
    const volScalarMomentFieldSet& moments,
    mappedPtrList<volScalarNode>& nodes
)
{
    // Recover reference to boundaryField of zero-order moment.
    // All moments will share the same BC types at a given boundary.
    const volScalarField::Boundary& bf = moments().boundaryField();

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            univariateMomentSet momentsToInvert
            (
                moments.size(),
                moments.support()
            );

            // Copying moments from a face
            forAll(momentsToInvert, momenti)
            {
                momentsToInvert[momenti]
                        = moments[momenti].boundaryField()[patchi][facei];
            }

            // Inverting moments for EQMOM
            momentInverter_->invert(momentsToInvert);

            // Recovering primary weights and abscissae from moment inverter
            const scalarList& pWeights(momentInverter_().primaryWeights());

            const scalarList& pAbscissae
            (
                momentInverter_().primaryAbscissae()
            );

            // Copying quadrature data to boundary face
            for (label pNodei = 0; pNodei < pWeights.size(); pNodei++)
            {
                volScalarNode& node = nodes[pNodei];

                node.primaryWeight().boundaryFieldRef()[patchi][facei]
                        = pWeights[pNodei];

                node.primaryAbscissae()[0].boundaryFieldRef()[patchi][facei]
                        = pAbscissae[pNodei];

                node.sigmas()[0].boundaryFieldRef()[patchi][facei]
                        = momentInverter_->sigma();

                for
                (
                    label sNodei = 0;
                    sNodei < node.nSecondaryNodes();
                    sNodei++
                )
                {
                    node.secondaryWeights()[0][sNodei].boundaryFieldRef()[patchi][facei]
                            = momentInverter_().secondaryWeights()[pNodei][sNodei];

                    node.secondaryAbscissae()[0][sNodei].boundaryFieldRef()[patchi][facei]
                            = momentInverter_().secondaryAbscissae()[pNodei][sNodei];
                }
            }
            for
            (
                label pNodei = pWeights.size();
                pNodei < nodes.size();
                pNodei++
            )
            {
                volScalarNode& node = nodes[pNodei];

                node.primaryWeight().boundaryFieldRef()[patchi][facei]
                        = 0.0;
                node.primaryAbscissae()[0].boundaryFieldRef()[patchi][facei]
                        = 0.0;

                node.sigmas()[0].boundaryFieldRef()[patchi][facei]
                        = 0.0;

                for
                (
                    label sNodei = 0;
                    sNodei < node.nSecondaryNodes();
                    sNodei++
                )
                {
                    node.secondaryWeights()[0][sNodei].boundaryFieldRef()[patchi][facei]
                            = 0.0;

                    node.secondaryAbscissae()[0][sNodei].boundaryFieldRef()[patchi][facei]
                            = 0.0;
                }
            }
        }
    }
}

bool Foam::extendedFieldMomentInversion::invertLocalMoments
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
        moments.support()
    );

    // Copying moment set from a cell to univariateMomentSet
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

    // Inverting moments and updating EQMOM
    momentInverter_().invert(momentsToInvert);

    // Recovering primary weights and abscissae from moment inverter
    const scalarList& pWeights(momentInverter_().primaryWeights());

    const scalarList& pAbscissae
    (
        momentInverter_().primaryAbscissae()
    );

    // Copying EQMOM quadrature to fields
    for (label pNodei = 0; pNodei < pWeights.size(); pNodei++)
    {
        volScalarNode& node(nodes[pNodei]);

        // Copy primary node
        node.primaryWeight()[celli] = pWeights[pNodei];
        node.primaryAbscissae()[0][celli] = pAbscissae[pNodei];

        // Copy secondary nodes
        PtrList<volScalarField>& sWeightFields(node.secondaryWeights()[0]);
        PtrList<volScalarField>& sAbscissaFields(node.secondaryAbscissae()[0]);

        const scalarRectangularMatrix& sWeights
        (
            momentInverter_().secondaryWeights()
        );

        const scalarRectangularMatrix& sAbscissae
        (
            momentInverter_().secondaryAbscissae()
        );

        for
        (
            label sNodei = 0;
            sNodei < nodes[0].nSecondaryNodes();
            sNodei++
        )
        {
            sWeightFields[sNodei][celli] = sWeights[pNodei][sNodei];
            sAbscissaFields[sNodei][celli] = sAbscissae[pNodei][sNodei];
        }

        // Copy sigma
        node.sigmas()[0][celli] = momentInverter_().sigma();
    }

    return true;
}

void Foam::extendedFieldMomentInversion::invert
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes
)
{
    NotImplemented;
}

void Foam::extendedFieldMomentInversion::invertBoundaryMoments
(
    const volVelocityMomentFieldSet& moments,
    mappedPtrList<volVelocityNode>& nodes
)
{
    NotImplemented;
}

bool Foam::extendedFieldMomentInversion::invertLocalMoments
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
