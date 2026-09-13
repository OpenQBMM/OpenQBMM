/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2026 Alberto Passalacqua
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

Application
    Test-TensorProductMomentInversion

Description
    Test the tensor product moment inversion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "mappedLists.H"
#include "supportType.H"
#include "TensorProductMomentInversion.H"
#include "Random.H"
#include "multivariateMomentTest.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Invert a moment set whose zero-order moment is below the smallest one
//  the inversion accepts.
//
//  Such a set carries no information, and the realizability check settles
//  it before any quadrature is built, so the inversion reports the failure
//  and leaves the quadrature at the zero it starts from. The caller is what
//  decides what to do with a cell like that, so it has to be told.
void testNegligibleMass
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
{
    Info<< "\n\nInverting a moment set of negligible mass" << endl;

    multivariateMomentSet moments
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(momentOrders[0].size(), supportType::R),
        SMALL,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        moments(momentOrders[mi]) = 0.1*SMALL;
    }

    multivariateMomentInversions::TensorProduct inverter
    (
        dict, momentOrders, nodeIndexes, velocityIndexes
    );

    if (inverter.invert(moments))
    {
        FatalErrorInFunction
            << "The inversion of a moment set of negligible mass was"
            << " reported as a success." << nl
            << "    Zero-order moment: " << moments(momentOrders[0]) << nl
            << "    Smallest accepted: " << inverter.smallM0() << nl
            << exit(FatalError);
    }

    Info<< "  the inversion reported the failure" << endl;

    forAll(inverter.weights(), nodei)
    {
        if (mag(inverter.weights()[nodei]) > SMALL)
        {
            FatalErrorInFunction
                << "A failed inversion left a weight behind." << nl
                << "    Node: " << nodei << nl
                << "    Weight: " << inverter.weights()[nodei] << nl
                << exit(FatalError);
        }
    }

    Info<< "  the quadrature it leaves behind is null" << endl;
}


//- Invert the moment set with the abscissae of one direction scaled, and
//  require the quadrature to follow: the weights and the abscissae of the
//  other directions unchanged, and those of the scaled direction scaled.
//
//  Each direction is inverted on its own, by its own quadrature, and its
//  zeta_k carry the units of its own abscissa. The inversion used to impose
//  on every direction the largest of their thresholds, floored at SMALL, so
//  a direction small enough in its own units, as a volume in cubic metres
//  is, was called degenerate however well spread it was.
//
//  knownFailure reports the outcome without asserting it, and warns when the
//  invariance starts holding, so that the pin is removed once it is fixed.
//  A failure may be a fatal error inside the inversion rather than a wrong
//  quadrature, so the inversion is run with exceptions on.
void testScaleInvariance
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes,
    const multivariateMomentSet& moments,
    const bool knownFailure = false
)
{
    const label scaledDim = 0;
    const scalar scale = 1.0e-18;

    Info<< "\n\nInverting the moment set with direction " << scaledDim
        << " scaled by " << scale << endl;

    multivariateMomentSet scaled
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(momentOrders[0].size(), supportType::R),
        SMALL,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scaled(momentOrder) =
            pow(scale, momentOrder[scaledDim])*moments(momentOrder);
    }

    const label nDims = momentOrders[0].size();
    const scalar tolerance = 1.0e-10;

    scalar worstWeight = VGREAT;
    vector worstAbscissa(VGREAT, VGREAT, VGREAT);
    string failure;

    const bool throwing = FatalError.throwing(true);

    try
    {
        multivariateMomentInversions::TensorProduct reference
        (
            dict, momentOrders, nodeIndexes, velocityIndexes
        );

        multivariateMomentInversions::TensorProduct inverter
        (
            dict, momentOrders, nodeIndexes, velocityIndexes
        );

        if (!reference.invert(moments) || !inverter.invert(scaled))
        {
            failure = "the inversion reported a failure";
        }
        else
        {
            // Each list is measured against the largest of its own entries,
            // so that a component that is zero is not asked to be reproduced
            // to no error
            scalar weightScale = 0;
            vector abscissaScale(Zero);

            forAll(nodeIndexes, nodei)
            {
                const labelList& nodeIndex = nodeIndexes[nodei];

                weightScale =
                    max(weightScale, mag(reference.weights()(nodeIndex)));

                for (label dimi = 0; dimi < nDims; dimi++)
                {
                    abscissaScale[dimi] =
                        max
                        (
                            abscissaScale[dimi],
                            mag(reference.velocityAbscissae()(nodeIndex)[dimi])
                        );
                }
            }

            worstWeight = 0;
            worstAbscissa = Zero;

            forAll(nodeIndexes, nodei)
            {
                const labelList& nodeIndex = nodeIndexes[nodei];

                worstWeight =
                    max
                    (
                        worstWeight,
                        mag
                        (
                            inverter.weights()(nodeIndex)
                          - reference.weights()(nodeIndex)
                        )/max(weightScale, SMALL)
                    );

                for (label dimi = 0; dimi < nDims; dimi++)
                {
                    const scalar factor = (dimi == scaledDim ? scale : 1.0);

                    worstAbscissa[dimi] =
                        max
                        (
                            worstAbscissa[dimi],
                            mag
                            (
                                inverter.velocityAbscissae()(nodeIndex)[dimi]
                               /factor
                              - reference.velocityAbscissae()(nodeIndex)[dimi]
                            )/max(abscissaScale[dimi], SMALL)
                        );
                }
            }
        }
    }
    catch (const Foam::error& err)
    {
        failure = err.message();
    }

    FatalError.throwing(throwing);

    const bool holds =
        failure.empty()
     && worstWeight <= tolerance
     && cmptMax(worstAbscissa) <= tolerance;

    if (failure.empty())
    {
        Info<< "  the weights differ by " << worstWeight
            << ", the abscissae of each direction by " << worstAbscissa
            << endl;
    }
    else
    {
        Info<< "  the inversion failed: " << failure.c_str() << endl;
    }

    if (knownFailure)
    {
        Info<< (holds ? "  (pinned)" : "  (known failure)") << endl;

        if (holds)
        {
            WarningInFunction
                << "Scaling one direction of the tensor product is pinned as a "
                << "known failure, and now holds." << nl
                << "    Remove the pin from the test, so that the defect is "
                << "caught if it returns." << nl << endl;
        }

        return;
    }

    if (!holds)
    {
        FatalErrorInFunction
            << "Scaling the abscissae of one direction did not scale the "
            << "quadrature with them." << nl
            << "    Direction scaled: " << scaledDim << ", by " << scale << nl
            << "    Largest difference in the weights: " << worstWeight << nl
            << "    Largest difference in the abscissae of each direction: "
            << worstAbscissa << nl
            << "    Tolerance: " << tolerance << nl
            << exit(FatalError);
    }

    Info<< "  the quadrature follows the scale" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "createFields.H"

    mappedList<scalarList> x
    (
        nNodes,
        nodeIndexes,
        scalarField(nDims, Zero)
    );

    mappedList<scalar> w(nNodes, nodeIndexes, 0.0);

    // A fixed seed, so that the moments the inversion is asked for are the
    // same on every run and on every machine
    Random rndGen(20260904);

    forAll(x, nodei)
    {
        w[nodei] = rndGen.sample01<scalar>();

        forAll(x[nodei], dimi)
        {
            x[nodei][dimi] = rndGen.sample01<scalar>();
        }
    }

    Info<< "Original moments:" << endl;

    multivariateMomentSet moments
    (
        nMoments,
        momentOrders,
        List<supportType>(momentOrders[0].size(), supportType::R),
        SMALL,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        moments(momentOrder) = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];
            scalar cmpt = w(nodeIndex);

            forAll(nodeIndex, dimi)
            {
                cmpt *= pow(x(nodeIndex)[dimi], momentOrder[dimi]);
            }

            moments(momentOrder) += cmpt;
        }

        Info<< "moment.";

        forAll(momentOrder, dimi)
        {
            Info<< momentOrder[dimi];
        }

        Info<< ": " << moments(momentOrder) << endl;
    }

    multivariateMomentInversions::TensorProduct momentInverter
    (
        quadratureProperties, momentOrders, nodeIndexes, velocityIndexes
    );

    Info<< "\nInverting moments" << endl;

    momentInverter.invert(moments);

    Info<< "\nReconstructed moments:" << endl;

    const mappedScalarList& weights = momentInverter.weights();
    const mappedList<scalarList>& abscissae = momentInverter.abscissae();

    const mappedVectorList& velocityAbscissae =
        momentInverter.velocityAbscissae();

    mappedList<scalar> newMoments(nMoments, momentOrders);

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        newMoments(momentOrder) = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            scalar cmpt = weights(nodeIndex);
            label vi = 0;
            label si = 0;

            for(label dimi = 0; dimi < momentOrder.size(); dimi++)
            {
                if (vi < velocityIndexes.size() && velocityIndexes[vi] == dimi)
                {
                     cmpt *=
                        pow
                        (
                            velocityAbscissae(nodeIndex)[vi],
                            momentOrder[dimi]
                        );
                    vi++;
                }
                else
                {
                    cmpt *= pow(abscissae(nodeIndex)[si], momentOrder[dimi]);
                    si++;
                }
            }

            newMoments(momentOrder) += cmpt;
        }

    }

    // The tensor product quadrature has a node for every combination of the
    // univariate abscissae, and solves for its weights from the mixed
    // moments, so every moment of the set is reproduced
    checkMomentConservation
    (
        newMoments,
        moments,
        momentOrders,
        1e-10,
        "TensorProduct"
    );


    testNegligibleMass
    (
        quadratureProperties,
        momentOrders,
        nodeIndexes,
        velocityIndexes
    );

    // Pinned: the weights are solved for from a tensor-product Vandermonde
    // matrix built from the abscissae in the units of the case, whose
    // entries span 1 to scale^2 for the scaled direction, and which the
    // elimination then finds singular. Solving it in the units of each
    // direction is what makes this hold, and the pin comes out with it.
    testScaleInvariance
    (
        quadratureProperties,
        momentOrders,
        nodeIndexes,
        velocityIndexes,
        moments,
        true                                // known failure, see above
    );

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
