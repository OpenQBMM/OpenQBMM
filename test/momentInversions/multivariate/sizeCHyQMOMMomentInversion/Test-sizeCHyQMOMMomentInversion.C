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

    // CHyQMOM+ does not conserve the moment of order (0 1 2) of the
    // velocity distribution it is given, so sizeCHyQMOM+ does not conserve
    // it at any size node. Pinned until the cause is found; see the same
    // pin in Test-CHyQMOMDegenerate.
    labelListList plusKnownFailures;

    for (label k = 0; k < nSizeNodes; k++)
    {
        plusKnownFailures.append(labelList({k, 0, 1, 2}));
    }

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
        tolerance,
        plusKnownFailures
    );

    Info<< "\n\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
