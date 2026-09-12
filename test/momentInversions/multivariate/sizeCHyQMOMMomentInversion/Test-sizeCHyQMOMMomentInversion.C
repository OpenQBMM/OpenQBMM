/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2026 Alberto Passalacqua
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
    Test-sizeCHyQMOMMomentInversion

Description
    Test the size conditioned velocity moment inversion classes.

    The distribution is built so that the method can represent it exactly:
    a size quadrature of as many nodes as the inversion builds, and, for
    each size node, a velocity quadrature of three nodes per direction,
    which is what CHyQMOM builds. Every moment the method controls then has
    to come back to round-off.

    A distribution drawn at random does not have that property. Its
    conditional velocity moments need not be realizable, and CHyQMOM
    corrects an unrealizable moment vector onto the boundary of the
    realizable region rather than reproducing it, so a third order moment
    comes back with a few percent of error and the method is behaving as
    designed. Such a test can only check that nothing crashed.

    The same distribution is then written at another scale, and the
    quadrature has to follow it. Scaling the size coordinate by s
    multiplies the moment of size order k by s^k, and has to leave the
    weights and the velocities where they are and multiply the size
    abscissae by s; scaling the measure by c multiplies every moment by c,
    and has to multiply the weights by c and leave everything else. Neither
    holds when a threshold inside the inversion is a fixed number compared
    with a quantity that carries the units of the abscissa or of the
    weight, which is what these two checks are for: the same population of
    droplets is the same population whether its size is written in metres,
    in cubic metres or in kilograms.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "IFstream.H"
#include "mappedLists.H"
#include "supportType.H"
#include "sizeCHyQMOMMomentInversions.H"
#include "multivariateMomentTest.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A distribution of a size and three velocity components, as a quadrature
//  of nSizeNodes size nodes, each carrying a velocity quadrature of three
//  nodes per direction
class sizeVelocityQuadrature
{
    //- Weights and abscissae of the size direction
    scalarList ws_;
    scalarList xs_;

    //- Weights of each velocity direction, which sum to one
    scalarList wv_;

    //- Abscissae of each velocity direction, of zero mean and unit variance
    scalarList xv_;


public:

    sizeVelocityQuadrature(const label nSizeNodes)
    :
        ws_(nSizeNodes, Zero),
        xs_(nSizeNodes, Zero),
        wv_({1.0/6.0, 2.0/3.0, 1.0/6.0}),
        xv_({-Foam::sqrt(scalar(3)), 0, Foam::sqrt(scalar(3))})
    {
        // Distinct, positive size abscissae, so that the size support is
        // R+ and the Vandermonde matrix the inversion builds from them is
        // not singular
        forAll(ws_, nodei)
        {
            ws_[nodei] = 0.5 + 0.25*nodei;
            xs_[nodei] = 0.4 + 0.7*nodei;
        }
    }

    //- Mean of velocity direction dimi at size node nodei. The velocity
    //  depends on the size, which is what the conditioning is for.
    scalar meanVelocity(const label nodei, const label dimi) const
    {
        static const scalarList slope({0.8, -0.5, 0.3});
        static const scalarList offset({-0.2, 0.6, 0.1});

        return offset[dimi] + slope[dimi]*xs_[nodei];
    }

    //- Standard deviation of velocity direction dimi at size node nodei
    scalar sigma(const label nodei, const label dimi) const
    {
        static const scalarList base({0.9, 1.3, 0.7});

        return base[dimi]*(1.0 + 0.2*nodei);
    }

    //- Moment of the given order
    scalar moment(const labelList& momentOrder, const label nVelocityDims)
    const
    {
        // Shear that correlates the velocity directions with one another
        // without changing the variance of any of them
        const scalar shearVU = 0.3;
        const scalar shearWU = -0.2;
        const scalar shearWV = 0.4;

        scalar m = 0.0;

        forAll(ws_, nodei)
        {
            scalar conditional = 0.0;

            forAll(xv_, i)
            {
                forAll(xv_, j)
                {
                    forAll(xv_, k)
                    {
                        const scalar u =
                            meanVelocity(nodei, 0) + sigma(nodei, 0)*xv_[i];

                        const scalar v =
                            meanVelocity(nodei, 1) + sigma(nodei, 1)*xv_[j]
                          + shearVU*u;

                        const scalar w =
                            meanVelocity(nodei, 2) + sigma(nodei, 2)*xv_[k]
                          + shearWU*u + shearWV*v;

                        const scalarList velocity({u, v, w});

                        scalar cmpt = wv_[i]*wv_[j]*wv_[k];

                        for (label dimi = 0; dimi < nVelocityDims; dimi++)
                        {
                            cmpt *=
                                pow(velocity[dimi], momentOrder[dimi + 1]);
                        }

                        conditional += cmpt;
                    }
                }
            }

            m += ws_[nodei]*pow(xs_[nodei], momentOrder[0])*conditional;
        }

        return m;
    }
};


//- The moments a sizeCHyQMOM inversion controls.
//
//  It builds nSizeNodes size nodes from the pure size moments, and then
//  recovers, for every moment order of the velocity method, one conditional
//  velocity moment per size node by inverting the Vandermonde matrix of the
//  size abscissae. So it consumes, and has to conserve, the pure size
//  moments and the mixed moments whose size order is below the number of
//  size nodes. Everything else in the moment set is closed rather than
//  conserved: a moment such as (2 0 1 2) has a velocity order that belongs
//  to CHyQMOM+ rather than to CHyQMOM, and is not asserted here.
template<class velocityInversion>
Foam::labelListList controlledMomentOrders
(
    const labelListList& momentOrders,
    const label nSizeMoments,
    const label nSizeNodes,
    const label nVelocityDims
)
{
    const labelListList velocityOrders
    (
        velocityInversion::getMomentOrders(nVelocityDims)
    );

    labelListList controlled;

    // The pure size moments
    for (label k = 0; k < nSizeMoments; k++)
    {
        labelList order(nVelocityDims + 1, 0);
        order[0] = k;
        controlled.append(order);
    }

    // The mixed moments of each size node
    for (label k = 0; k < nSizeNodes; k++)
    {
        forAll(velocityOrders, vi)
        {
            labelList order(nVelocityDims + 1, 0);
            order[0] = k;

            bool pureSize = true;

            for (label dimi = 0; dimi < nVelocityDims; dimi++)
            {
                order[dimi + 1] = velocityOrders[vi][dimi];

                if (order[dimi + 1] != 0)
                {
                    pureSize = false;
                }
            }

            // Only assert the moments the set actually carries, and leave
            // the pure size moments to the list above
            bool present = false;

            forAll(momentOrders, mi)
            {
                if (momentOrders[mi] == order)
                {
                    present = true;
                    break;
                }
            }

            if (!pureSize && present)
            {
                controlled.append(order);
            }
        }
    }

    return controlled;
}


//- Invert the moments and assert that the method conserves its own
template<class inversionType, class velocityInversion>
void testInversion
(
    const word& what,
    const dictionary& dict,
    const multivariateMomentSet& moments,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes,
    const label nDims,
    const label nSizeMoments,
    const label nSizeNodes,
    const scalar tolerance,
    const labelListList& knownFailures = labelListList()
)
{
    Info<< "\n\nInverting moments with " << what << endl;

    inversionType inverter(dict, momentOrders, nodeIndexes, velocityIndexes);

    if (!inverter.invert(moments))
    {
        FatalErrorInFunction
            << "The inversion with " << what << " failed." << nl
            << exit(FatalError);
    }

    const mappedScalarList& weights = inverter.weights();
    const mappedList<scalarList>& sizeAbscissae = inverter.abscissae();
    const mappedVectorList& velocityAbscissae = inverter.velocityAbscissae();

    mappedList<scalar> reconstructed(momentOrders.size(), momentOrders, Zero);

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scalar m = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            scalar cmpt = weights(nodeIndex);

            cmpt *= pow(sizeAbscissae(nodeIndex)[0], momentOrder[0]);

            for (label dimi = 0; dimi < nDims - 1; dimi++)
            {
                cmpt *=
                    pow
                    (
                        velocityAbscissae(nodeIndex)[dimi],
                        momentOrder[dimi + 1]
                    );
            }

            m += cmpt;
        }

        reconstructed(momentOrder) = m;
    }

    checkMomentConservation
    (
        reconstructed,
        moments,
        controlledMomentOrders<velocityInversion>
        (
            momentOrders,
            nSizeMoments,
            nSizeNodes,
            nDims - 1
        ),
        tolerance,
        what,
        knownFailures
    );
}


//- What an inversion returns, copied out of the inverter so that two of
//  them can be compared
struct quadratureOf
{
    bool inverted;
    scalarList w;
    scalarList x;
    List<vector> U;
};


//- Invert a moment set and copy the quadrature out
template<class inversionType>
quadratureOf invertOnce
(
    const dictionary& dict,
    const multivariateMomentSet& moments,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
{
    inversionType inverter(dict, momentOrders, nodeIndexes, velocityIndexes);

    quadratureOf q;

    q.inverted = inverter.invert(moments);
    q.w.setSize(nodeIndexes.size(), Zero);
    q.x.setSize(nodeIndexes.size(), Zero);
    q.U.setSize(nodeIndexes.size(), Zero);

    if (!q.inverted)
    {
        return q;
    }

    forAll(nodeIndexes, nodei)
    {
        const labelList& nodeIndex = nodeIndexes[nodei];

        q.w[nodei] = inverter.weights()(nodeIndex);
        q.x[nodei] = inverter.abscissae()(nodeIndex)[0];
        q.U[nodei] = inverter.velocityAbscissae()(nodeIndex);
    }

    return q;
}


//- Compare two quadratures, each list against the largest of its own
//  entries, so that a component that is zero is not asked to be reproduced
//  to no error at all.
//
//  knownFailure reports the outcome without asserting it, for an invariance
//  that is known not to hold yet. One that starts holding is reported too,
//  because the pin has to be removed once it is fixed: left in place it
//  would hide the defect coming back.
bool compareQuadratures
(
    const quadratureOf& computed,
    const quadratureOf& expected,
    const scalar tolerance,
    const word& what,
    const bool knownFailure = false
)
{
    Info<< "\nVerifying that " << what << endl;

    bool holds = computed.inverted && expected.inverted;

    if (!holds)
    {
        Info<< "  the inversion failed" << endl;
    }
    else
    {
        scalar wScale = 0.0;
        scalar xScale = 0.0;
        scalar UScale = 0.0;

        forAll(expected.w, nodei)
        {
            wScale = max(wScale, mag(expected.w[nodei]));
            xScale = max(xScale, mag(expected.x[nodei]));
            UScale = max(UScale, cmptMax(cmptMag(expected.U[nodei])));
        }

        const scalar wBound = tolerance*max(wScale, SMALL);
        const scalar xBound = tolerance*max(xScale, SMALL);
        const scalar UBound = tolerance*max(UScale, SMALL);

        scalar wWorst = 0.0;
        scalar xWorst = 0.0;
        scalar UWorst = 0.0;

        forAll(expected.w, nodei)
        {
            wWorst = max(wWorst, mag(computed.w[nodei] - expected.w[nodei]));
            xWorst = max(xWorst, mag(computed.x[nodei] - expected.x[nodei]));

            UWorst =
                max
                (
                    UWorst,
                    cmptMax(cmptMag(computed.U[nodei] - expected.U[nodei]))
                );
        }

        Info<< "  weights differ by at most " << wWorst
            << ", of " << wScale << nl
            << "  sizes differ by at most " << xWorst
            << ", of " << xScale << nl
            << "  velocities differ by at most " << UWorst
            << ", of " << UScale << endl;

        holds = (wWorst <= wBound && xWorst <= xBound && UWorst <= UBound);
    }

    if (knownFailure)
    {
        Info<< (holds ? "  (pinned)" : "  (known failure)") << endl;

        if (holds)
        {
            WarningInFunction
                << what << " is pinned as a known failure, and now holds."
                << nl
                << "    Remove the pin from the test, so that the defect is"
                << " caught if it returns." << nl << endl;
        }

        return true;
    }

    if (!holds)
    {
        FatalErrorInFunction
            << "It is not the case that " << what << "." << nl
            << exit(FatalError);
    }

    Info<< "  it does" << endl;

    return true;
}


//- The same moment set with the size coordinate scaled by s and the whole
//  measure by c
multivariateMomentSet scaledMoments
(
    const multivariateMomentSet& moments,
    const labelListList& momentOrders,
    const label nDims,
    const scalar s,
    const scalar c,
    const scalar smallM0
)
{
    multivariateMomentSet scaled
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(nDims, supportType::R),
        smallM0,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scaled(momentOrder) =
            c*pow(s, momentOrder[0])*moments(momentOrder);
    }

    return scaled;
}


int main()
{
    #include "createFields.H"

    // The number of pure size moments the set carries, and the number of
    // size nodes the inversion builds from them
    label nSizeMoments = 0;

    forAll(momentOrders, mi)
    {
        nSizeMoments = max(nSizeMoments, momentOrders[mi][0] + 1);
    }

    label nSizeNodes = 0;

    forAll(nodeIndexes, nodei)
    {
        nSizeNodes = max(nSizeNodes, nodeIndexes[nodei][0] + 1);
    }

    Info<< "Size moments: " << nSizeMoments
        << ", size nodes: " << nSizeNodes << nl << endl;

    const sizeVelocityQuadrature source(nSizeNodes);

    multivariateMomentSet moments
    (
        nMoments,
        momentOrders,
        List<supportType>(nDims, supportType::R),
        SMALL,
        SMALL
    );

    Info<< "Moments of the quadrature:" << endl;

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        moments(momentOrder) = source.moment(momentOrder, nDims - 1);

        Info<< "  " << momentName(momentOrder) << ": "
            << moments(momentOrder) << endl;
    }

    const scalar tolerance = 1e-10;

    testInversion
    <
        multivariateMomentInversions::sizeCHyQMOM,
        multivariateMomentInversions::CHyQMOM
    >
    (
        "sizeCHyQMOM",
        quadratureProperties,
        moments,
        momentOrders,
        nodeIndexes,
        velocityIndexes,
        nDims,
        nSizeMoments,
        nSizeNodes,
        tolerance
    );

    testInversion
    <
        multivariateMomentInversions::sizeCHyQMOMPlus,
        multivariateMomentInversions::CHyQMOMPlus
    >
    (
        "sizeCHyQMOMPlus",
        quadratureProperties,
        moments,
        momentOrders,
        nodeIndexes,
        velocityIndexes,
        nDims,
        nSizeMoments,
        nSizeNodes,
        tolerance
    );


    // * * * * * * * * * * * * * * The two scales * * * * * * * * * * * * //

    // The measure written small. A weight is then a small number in its own
    // units, as a volume fraction of particles of a few nanometres is, and
    // the threshold the moment of order zero is cut at has to be lowered
    // with it: that entry is written here, where the inversion reads it.
    {
        const scalar c = 1.0e-20;
        const scalar smallM0 = 1.0e-40;

        dictionary smallMeasure(quadratureProperties);
        smallMeasure.set("smallM0", smallM0);

        const quadratureOf reference
        (
            invertOnce<multivariateMomentInversions::sizeCHyQMOM>
            (
                smallMeasure, moments, momentOrders, nodeIndexes,
                velocityIndexes
            )
        );

        quadratureOf expected(reference);

        forAll(expected.w, nodei)
        {
            expected.w[nodei] *= c;
        }

        const quadratureOf computed
        (
            invertOnce<multivariateMomentInversions::sizeCHyQMOM>
            (
                smallMeasure,
                scaledMoments
                (
                    moments, momentOrders, nDims, 1.0, c, smallM0
                ),
                momentOrders, nodeIndexes, velocityIndexes
            )
        );

        compareQuadratures
        (
            computed,
            expected,
            tolerance,
            "scaling the measure by " + Foam::name(c)
          + " scales the weights by it and moves nothing else"
        );
    }

    // The size coordinate written in smaller units, as a volume in cubic
    // metres is against a diameter in metres. What makes this hold is that
    // the size moments are normalised by their own mean before they are
    // inverted: the zeta_k the realizability check compares with smallZeta
    // carry the units of the abscissa, so without that a size small enough
    // in them is declared degenerate however well spread it is.
    {
        const scalar sizeScale = 1.0e-18;

        const quadratureOf reference
        (
            invertOnce<multivariateMomentInversions::sizeCHyQMOM>
            (
                quadratureProperties, moments, momentOrders, nodeIndexes,
                velocityIndexes
            )
        );

        quadratureOf expected(reference);

        forAll(expected.x, nodei)
        {
            expected.x[nodei] *= sizeScale;
        }

        const quadratureOf computed
        (
            invertOnce<multivariateMomentInversions::sizeCHyQMOM>
            (
                quadratureProperties,
                scaledMoments
                (
                    moments, momentOrders, nDims, sizeScale, 1.0, SMALL
                ),
                momentOrders, nodeIndexes, velocityIndexes
            )
        );

        compareQuadratures
        (
            computed,
            expected,
            tolerance,
            "scaling the size coordinate by " + Foam::name(sizeScale)
          + " scales the size abscissae by it and moves nothing else"
        );
    }

    Info<< "\n\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
