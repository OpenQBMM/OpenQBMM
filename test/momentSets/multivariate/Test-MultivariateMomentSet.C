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
    Test-MultivariateMomentSet

Description
    Test the multivariateMomentSet class.

    What the class carries is a set of moment values addressed by their
    order, so what the tests turn on is the addressing: a moment written
    with an order has to be the one read back with that order, whatever
    position it was declared at. The orders are packed into a single label
    key, and the packing is what the multivariate quadrature rests on.

    The rest is what the class refuses. Its constructor validates, so an
    order of the wrong number of components, a moment that is not finite,
    and two orders that land on the same key are all refused there rather
    than carried into an inversion.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "labelList.H"
#include "scalarList.H"
#include "supportType.H"
#include "multivariateMomentSet.H"
#include "CHyQMOMMomentInversion.H"
#include "CHyQMOMPlusMomentInversion.H"

#include <limits>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nTested = 0;

void check(const bool condition, const string& what)
{
    nTested++;

    if (!condition)
    {
        FatalErrorInFunction
            << "Failed: " << what << nl
            << exit(FatalError);
    }

    Info<< "  OK: " << what << endl;
}


void checkEqual
(
    const scalar computed,
    const scalar expected,
    const string& what
)
{
    check(mag(computed - expected) <= SMALL*max(mag(expected), 1.0), what);
}


//- Run f and require that it is refused with a fatal error. The moment set
//  validates in its constructor, so this is how the refusals are reached
//  without ending the test.
template<class Function>
void checkRefused(Function f, const string& what)
{
    nTested++;

    const bool oldThrowing = FatalError.throwing(true);
    bool refused = false;

    try
    {
        f();
    }
    catch (const Foam::error&)
    {
        refused = true;
    }

    FatalError.throwing(oldThrowing);

    if (!refused)
    {
        FatalErrorInFunction
            << "Accepted what it should refuse: " << what << nl
            << exit(FatalError);
    }

    Info<< "  OK: refused " << what << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Orders of a three-dimensional set, deliberately not in a sorted order,
//  so that a test that passes by reading the position rather than the key
//  is caught
labelListList threeDimOrders()
{
    return labelListList
    ({
        {0, 0, 0},
        {1, 1, 0},
        {0, 1, 0},
        {2, 0, 0},
        {0, 0, 1},
        {1, 0, 1},
        {0, 2, 0},
        {1, 0, 0},
        {0, 0, 2},
        {0, 1, 1}
    });
}


List<supportType> threeDimSupports()
{
    return List<supportType>(3, supportType::R);
}


//- A value that depends on the order, so that a moment read with the wrong
//  key is a different number
scalar valueOf(const labelList& order)
{
    return 1.0 + order[0] + 10.0*order[1] + 100.0*order[2];
}


void testConstructionFromSize()
{
    Info<< "\nConstructing from a size and an initial value" << endl;

    const labelListList orders(threeDimOrders());

    multivariateMomentSet moments
    (
        orders.size(), orders, threeDimSupports(), SMALL, SMALL, 3.5
    );

    check(moments.nMoments() == orders.size(), "the number of moments");
    check(moments.nDimensions() == 3, "the number of dimensions");
    check(moments.momentOrders() == orders, "the moment orders");
    check(moments.supports().size() == 3, "the number of supports");

    bool allSet = true;

    forAll(moments, mi)
    {
        allSet = allSet && (mag(moments[mi] - 3.5) <= SMALL);
    }

    check(allSet, "every moment carries the initial value");
}


void testConstructionFromList()
{
    Info<< "\nConstructing from a list of values" << endl;

    const labelListList orders(threeDimOrders());

    scalarList values(orders.size(), Zero);

    forAll(orders, mi)
    {
        values[mi] = valueOf(orders[mi]);
    }

    multivariateMomentSet moments
    (
        values, orders, threeDimSupports(), SMALL, SMALL
    );

    bool inOrder = true;

    forAll(orders, mi)
    {
        inOrder = inOrder && (mag(moments[mi] - values[mi]) <= SMALL);
    }

    check(inOrder, "the values are kept in the order they were given");

    bool byKey = true;

    forAll(orders, mi)
    {
        byKey =
            byKey && (mag(moments(orders[mi]) - valueOf(orders[mi])) <= SMALL);
    }

    check(byKey, "every moment is read back with its own order");
}


//- The heart of it: the same orders declared in two different sequences
//  have to address the same moments
void testAddressingIsIndependentOfDeclarationOrder()
{
    Info<< "\nAddressing a moment by its order rather than its position"
        << endl;

    const labelListList orders(threeDimOrders());

    labelListList shuffled(orders.size());

    forAll(orders, mi)
    {
        shuffled[mi] = orders[orders.size() - 1 - mi];
    }

    multivariateMomentSet a
    (
        orders.size(), orders, threeDimSupports(), SMALL, SMALL
    );
    multivariateMomentSet b
    (
        shuffled.size(), shuffled, threeDimSupports(), SMALL, SMALL
    );

    // Written by order into both
    forAll(orders, mi)
    {
        a(orders[mi]) = valueOf(orders[mi]);
        b(orders[mi]) = valueOf(orders[mi]);
    }

    bool agree = true;

    forAll(orders, mi)
    {
        agree = agree && (mag(a(orders[mi]) - b(orders[mi])) <= SMALL);
    }

    check(agree, "the two sets agree when read by order");

    // and the positions really do differ, so the test above has teeth
    check
    (
        mag(a[0] - b[0]) > SMALL,
        "the two sets differ when read by position"
    );
}


void testVariadicAccess()
{
    Info<< "\nAddressing with the orders given one by one" << endl;

    const labelListList orders(threeDimOrders());

    multivariateMomentSet moments
    (
        orders.size(), orders, threeDimSupports(), SMALL, SMALL
    );

    forAll(orders, mi)
    {
        moments(orders[mi]) = valueOf(orders[mi]);
    }

    checkEqual(moments(1, 1, 0), valueOf({1, 1, 0}), "moment (1 1 0)");
    checkEqual(moments(0, 0, 2), valueOf({0, 0, 2}), "moment (0 0 2)");
    checkEqual(moments(2, 0, 0), valueOf({2, 0, 0}), "moment (2 0 0)");
    checkEqual(moments(0, 1, 1), valueOf({0, 1, 1}), "moment (0 1 1)");
}


void testMomentMap()
{
    Info<< "\nThe map from an order to the position it was declared at"
        << endl;

    const labelListList orders(threeDimOrders());

    multivariateMomentSet moments
    (
        orders.size(), orders, threeDimSupports(), SMALL, SMALL
    );

    check
    (
        moments.momentMap().size() == orders.size(),
        "the map holds one key for every moment"
    );

    bool found = true;

    forAll(orders, mi)
    {
        found = found && moments.found(orders[mi]);
    }

    check(found, "every declared order is found");

    check
    (
        !moments.found(labelList({3, 0, 0})),
        "an order that was not declared is not found"
    );
}


void testRefusals()
{
    Info<< "\nWhat the moment set refuses" << endl;

    const labelListList orders(threeDimOrders());

    checkRefused
    (
        []()
        {
            labelListList bad(threeDimOrders());
            bad[2] = labelList({0, 1});          // two components, not three

            multivariateMomentSet moments
            (
                bad.size(), bad, threeDimSupports(), SMALL, SMALL
            );
        },
        "an order with the wrong number of components"
    );

    checkRefused
    (
        [&orders]()
        {
            scalarList values(orders.size(), Zero);
            values[3] = std::numeric_limits<scalar>::quiet_NaN();

            multivariateMomentSet moments
            (
                values, orders, threeDimSupports(), SMALL, SMALL
            );
        },
        "a moment that is not a number"
    );

    checkRefused
    (
        [&orders]()
        {
            scalarList values(orders.size(), Zero);
            values[1] = std::numeric_limits<scalar>::infinity();

            multivariateMomentSet moments
            (
                values, orders, threeDimSupports(), SMALL, SMALL
            );
        },
        "a moment that is not finite"
    );

    checkRefused
    (
        []()
        {
            labelListList repeated(threeDimOrders());
            repeated[4] = repeated[1];           // the same order twice

            multivariateMomentSet moments
            (
                repeated.size(), repeated, threeDimSupports(), SMALL, SMALL
            );
        },
        "the same order declared twice"
    );

    // The orders are packed into one label, a digit to a dimension, so an
    // order of ten or more in a dimension lands on the key of another. The
    // set refuses that rather than carrying two moments that address the
    // same value.
    checkRefused
    (
        []()
        {
            labelListList colliding
            ({
                {0, 0}, {1, 0}, {0, 10}
            });

            multivariateMomentSet moments
            (
                colliding.size(),
                colliding,
                List<supportType>(2, supportType::R),
                SMALL,
                SMALL
            );
        },
        "two orders that pack to the same key"
    );
}


void testSetSize()
{
    Info<< "\nResizing" << endl;

    const labelListList orders(threeDimOrders());

    multivariateMomentSet moments
    (
        orders.size(), orders, threeDimSupports(), SMALL, SMALL
    );

    labelListList smaller
    ({
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}
    });

    moments.setSize(smaller.size(), smaller);

    check(moments.nMoments() == smaller.size(), "the new number of moments");
    check(moments.momentOrders() == smaller, "the new moment orders");

    forAll(smaller, mi)
    {
        moments(smaller[mi]) = valueOf(smaller[mi]);
    }

    bool byKey = true;

    forAll(smaller, mi)
    {
        byKey =
            byKey
         && (mag(moments(smaller[mi]) - valueOf(smaller[mi])) <= SMALL);
    }

    check(byKey, "the moments are addressed by their new orders");

    checkRefused
    (
        [&moments]()
        {
            labelListList tooFew({{0, 0, 0}});
            multivariateMomentSet copy(moments);
            copy.setSize(1, tooFew);
        },
        "a resize to fewer than two moments"
    );
}


//- The inversions describe their own moment sets for one, two and three
//  dimensions. sizeCHyQMOM hands whatever the case declares straight to
//  them, so a distribution of four velocity dimensions - a dictionary that
//  gives four abscissae the dimensions of a velocity, which is a mistake
//  easily made with a size direction - used to be answered with an empty
//  set of orders rather than refused.
void testUnsupportedDimensions()
{
    Info<< "\nWhat the inversions refuse to describe" << endl;

    checkRefused
    (
        []() { multivariateMomentInversions::CHyQMOM::getNMoments(4); },
        "CHyQMOM asked for the number of moments of four dimensions"
    );

    checkRefused
    (
        []() { multivariateMomentInversions::CHyQMOM::getMomentOrders(4); },
        "CHyQMOM asked for the moment orders of four dimensions"
    );

    checkRefused
    (
        []() { multivariateMomentInversions::CHyQMOM::getNNodes(0); },
        "CHyQMOM asked for the number of nodes of no dimension"
    );

    checkRefused
    (
        []() { multivariateMomentInversions::CHyQMOM::getNodeIndexes(4); },
        "CHyQMOM asked for the node indexes of four dimensions"
    );

    checkRefused
    (
        []() { multivariateMomentInversions::CHyQMOMPlus::getNMoments(4); },
        "CHyQMOMPlus asked for the number of moments of four dimensions"
    );

    checkRefused
    (
        []() { multivariateMomentInversions::CHyQMOMPlus::getMomentOrders(4); },
        "CHyQMOMPlus asked for the moment orders of four dimensions"
    );

    checkRefused
    (
        []() { multivariateMomentInversions::CHyQMOMPlus::getNNodes(4); },
        "CHyQMOMPlus asked for the number of nodes of four dimensions"
    );

    checkRefused
    (
        []() { multivariateMomentInversions::CHyQMOMPlus::getNodeIndexes(4); },
        "CHyQMOMPlus asked for the node indexes of four dimensions"
    );
}


int main()
{
    testConstructionFromSize();
    testConstructionFromList();
    testAddressingIsIndependentOfDeclarationOrder();
    testVariadicAccess();
    testMomentMap();
    testRefusals();
    testSetSize();
    testUnsupportedDimensions();

    Info<< "\n" << nTested << " checks passed.\n" << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
