/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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
    Test-CHyQMOMDegenerate

Description
    Test the three-dimensional CHyQMOM and CHyQMOM+ inversions on velocity
    distributions with a degenerate direction.

    invert3D takes a different branch for each combination of directions
    whose variance falls below varMin, and each branch builds the velocity
    abscissae of the degenerate directions on its own. A direction with no
    variance still has a mean, and the quadrature has to carry it: the
    moment of first order in that direction is not zero.

    The distributions here are built from a quadrature of three nodes per
    active direction and a single node per degenerate one, so that the
    moments of the method are the moments of a distribution it can
    represent, and every one of them has to come back. A degenerate
    direction is given a mean far from zero, which is what tells a branch
    that carries the mean apart from one that leaves the abscissa at the
    value reset() wrote.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "IOstreams.H"
#include "mappedLists.H"
#include "supportType.H"
#include "multivariateMomentSet.H"
#include "CHyQMOMMomentInversion.H"
#include "CHyQMOMPlusMomentInversion.H"
#include "multivariateMomentTest.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A quadrature of the velocity distribution, as the tensor product of the
//  three directions. A direction of a single node is degenerate: it carries
//  a mean and no variance.
class velocityQuadrature
{
public:

    //- Weights of each direction, which sum to one
    scalarListList w;

    //- Abscissae of each direction
    scalarListList x;

    //- Shear of the second and third direction onto the first, and of the
    //  third onto the second, which correlates the directions without
    //  changing the variance of any of them
    scalar shearVU;
    scalar shearWU;
    scalar shearWV;

    //- Number density
    scalar m0;

    velocityQuadrature
    (
        const scalarListList& weights,
        const scalarListList& abscissae,
        const scalar shearVUIn = 0,
        const scalar shearWUIn = 0,
        const scalar shearWVIn = 0,
        const scalar m0In = 1.5
    )
    :
        w(weights),
        x(abscissae),
        shearVU(shearVUIn),
        shearWU(shearWUIn),
        shearWV(shearWVIn),
        m0(m0In)
    {}

    //- Moment of the given order
    scalar moment(const labelList& momentOrder) const
    {
        scalar m = 0.0;

        forAll(x[0], i)
        {
            forAll(x[1], j)
            {
                forAll(x[2], k)
                {
                    const scalar u = x[0][i];
                    const scalar v = x[1][j] + shearVU*u;
                    const scalar wc = x[2][k] + shearWU*u + shearWV*v;

                    m += m0*w[0][i]*w[1][j]*w[2][k]
                        *pow(u, momentOrder[0])
                        *pow(v, momentOrder[1])
                        *pow(wc, momentOrder[2]);
                }
            }
        }

        return m;
    }
};


//- Invert the moments of the quadrature and assert that every moment of the
//  method comes back
template<class inversionType>
void testInversion
(
    const word& what,
    const velocityQuadrature& source,
    const scalar tolerance,
    const labelListList& knownFailures = labelListList()
)
{
    const labelListList momentOrders(inversionType::getMomentOrders(3));
    const labelListList nodeIndexes(inversionType::getNodeIndexes(3));

    multivariateMomentSet moments
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(3, supportType::R),
        SMALL,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        moments(momentOrders[mi]) = source.moment(momentOrders[mi]);
    }

    dictionary dict;

    inversionType inverter(dict, momentOrders, nodeIndexes, {0, 1, 2});

    if (!inverter.invert(moments))
    {
        FatalErrorInFunction
            << "The inversion of " << what << " failed." << nl
            << exit(FatalError);
    }

    // Rebuild the moments from the quadrature the inversion returned
    multivariateMomentSet reconstructed
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(3, supportType::R),
        SMALL,
        SMALL
    );

    const mappedScalarList& weights = inverter.weights();
    const mappedVectorList& abscissae = inverter.velocityAbscissae();

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scalar m = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            scalar cmpt = weights(nodeIndex);

            for (label dimi = 0; dimi < 3; dimi++)
            {
                cmpt *= pow(abscissae(nodeIndex)[dimi], momentOrder[dimi]);
            }

            m += cmpt;
        }

        reconstructed(momentOrder) = m;
    }

    checkMomentConservation
    (
        reconstructed,
        moments,
        momentOrders,
        tolerance,
        what,
        knownFailures
    );
}


//- Run a case through both inversions
void testCase
(
    const word& name,
    const velocityQuadrature& source,
    const scalar tolerance,
    const labelListList& plusKnownFailures = labelListList()
)
{
    Info<< "\n\nTesting " << name << endl;

    testInversion<multivariateMomentInversions::CHyQMOM>
    (
        "CHyQMOM, " + name,
        source,
        tolerance
    );

    testInversion<multivariateMomentInversions::CHyQMOMPlus>
    (
        "CHyQMOMPlus, " + name,
        source,
        tolerance,
        plusKnownFailures
    );
}


//- Invert a two-dimensional velocity distribution whose first direction
//  carries no variance.
//
//  invert3D never hands invert2D a degenerate first direction: it settles
//  each degenerate combination in a branch of its own. The branch of
//  invert2D that takes one, the one that writes wDir1 directly, is
//  therefore reached only when the distribution itself has two velocity
//  dimensions, as the rotating cylinder and jet crossing cases have. It
//  had no coverage before this.
//
//  That branch writes wDir1[1] into a list it sized from what the
//  univariate inversion returned, and CHyQMOM reads three nodes from that
//  inversion everywhere without checking. It is safe only because CHyQMOM
//  always uses hyperbolicMomentInversion, which fixes nNodes to three; the
//  base class contract allows one node for a degenerate set and none for a
//  vanishing one. This case drives the degenerate input through it, so that
//  the day the univariate inversion is made selectable the assumption is
//  not silently broken.
template<class inversionType>
void testTwoDimensionalDegenerate(const word& what)
{
    Info<< "\n\nTesting " << what
        << ", two dimensions with a degenerate first direction" << endl;

    const labelListList momentOrders(inversionType::getMomentOrders(2));
    const labelListList nodeIndexes(inversionType::getNodeIndexes(2));

    multivariateMomentSet moments
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(2, supportType::R),
        SMALL,
        SMALL
    );

    // A delta at the origin in the first direction, and three nodes of
    // unit variance in the second. Every moment of the first direction is
    // exactly zero, so its central variance is exactly zero rather than the
    // round-off of a difference.
    moments(0, 0) = 1.0;
    moments(1, 0) = 0.0;
    moments(2, 0) = 0.0;
    moments(3, 0) = 0.0;
    moments(4, 0) = 0.0;
    moments(1, 1) = 0.0;
    moments(0, 1) = 0.0;
    moments(0, 2) = 1.0;
    moments(0, 3) = 0.0;
    moments(0, 4) = 3.0;

    dictionary dict;
    inversionType inverter(dict, momentOrders, nodeIndexes, {0, 1});

    inverter.invert(moments);

    const mappedScalarList& weights = inverter.weights();
    const mappedVectorList& abscissae = inverter.velocityAbscissae();

    mappedList<scalar> reconstructed(momentOrders.size(), momentOrders, Zero);

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        scalar m = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];
            scalar cmpt = weights(nodeIndex);

            for (label dimi = 0; dimi < 2; dimi++)
            {
                cmpt *= pow(abscissae(nodeIndex)[dimi], momentOrder[dimi]);
            }

            m += cmpt;
        }

        reconstructed(momentOrder) = m;
    }

    checkMomentConservation
    (
        reconstructed,
        moments,
        momentOrders,
        1e-10,
        what + ", two dimensions with a degenerate first direction"
    );
}


int main()
{
    // Three nodes of unit variance and no skewness, so that a direction
    // built from them is realizable up to fourth order
    const scalar root3 = Foam::sqrt(scalar(3));
    const scalarList wActive({1.0/6.0, 2.0/3.0, 1.0/6.0});
    const scalarList xActive({-root3, scalar(0), root3});

    // A degenerate direction is a single node away from the origin, so that
    // an abscissa left at zero is not mistaken for the mean it should carry
    const scalarList wFrozen({1.0});

    const scalar tolerance = 1e-10;

    // Every direction active, which is the branch the other tests cover.
    //
    // CHyQMOM+ does not conserve the moment of order (0 1 2) here. The
    // distribution is symmetric in every direction, so that moment is zero,
    // and the inversion returns minus the moment of order (1 1 0). It is
    // pinned rather than asserted until the cause is found: it is a defect
    // of CHyQMOM+ in its own right and not a property of the degenerate
    // directions this test covers. CHyQMOM conserves it.
    testCase
    (
        "no degenerate direction",
        velocityQuadrature
        (
            {wActive, wActive, wActive},
            {xActive, xActive, xActive},
            0.3, -0.2, 0.4
        ),
        tolerance,
        {{0, 1, 2}}
    );

    // One degenerate direction at a time
    testCase
    (
        "degenerate z direction",
        velocityQuadrature
        (
            {wActive, wActive, wFrozen},
            {xActive, xActive, {2.5}},
            0.3
        ),
        tolerance
    );

    testCase
    (
        "degenerate x direction",
        velocityQuadrature
        (
            {wFrozen, wActive, wActive},
            {{-1.75}, xActive, xActive},
            0, 0, 0.4
        ),
        tolerance
    );

    testCase
    (
        "degenerate y direction",
        velocityQuadrature
        (
            {wActive, wFrozen, wActive},
            {xActive, {3.25}, xActive},
            0, -0.2
        ),
        tolerance
    );

    // Two degenerate directions
    testCase
    (
        "degenerate x and y directions",
        velocityQuadrature
        (
            {wFrozen, wFrozen, wActive},
            {{-1.75}, {3.25}, xActive}
        ),
        tolerance
    );

    // Two degenerate directions whose single node sits *at the origin*.
    //
    // The cases above put the frozen node away from zero, which leaves the
    // central variance of that direction at the round-off of mu*mu - mu*mu,
    // and so on either side of zero. Here every moment of the frozen
    // directions is exactly zero, so the central variance is exactly zero
    // and the realizability correction zeroes the third and fourth central
    // moments with it. That is the input that reports the univariate set
    // degenerate, which is the one worth pinning.
    testCase
    (
        "degenerate x and y directions at the origin",
        velocityQuadrature
        (
            {wFrozen, wFrozen, wActive},
            {{0.0}, {0.0}, xActive}
        ),
        tolerance
    );

    testCase
    (
        "degenerate x and z directions",
        velocityQuadrature
        (
            {wFrozen, wActive, wFrozen},
            {{-1.75}, xActive, {2.5}}
        ),
        tolerance
    );

    // Every direction degenerate, which leaves a single velocity
    testCase
    (
        "every direction degenerate",
        velocityQuadrature
        (
            {wFrozen, wFrozen, wFrozen},
            {{-1.75}, {3.25}, {2.5}}
        ),
        tolerance
    );

    testTwoDimensionalDegenerate<multivariateMomentInversions::CHyQMOM>
    (
        "CHyQMOM"
    );

    testTwoDimensionalDegenerate<multivariateMomentInversions::CHyQMOMPlus>
    (
        "CHyQMOMPlus"
    );

    Info<< "\n\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
