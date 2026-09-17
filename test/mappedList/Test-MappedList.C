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
    Test-MappedList

Description
    Test the mappedList class.

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "mappedList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkScalar
(
    const word& description,
    const scalar computed,
    const scalar expected
)
{
    if (mag(computed - expected) > SMALL)
    {
        FatalErrorInFunction
            << description << " does not match: " << endl
            << "  Expected value: " << expected << endl
            << "  Computed value: " << computed << endl
            << exit(FatalError);
    }
}

void checkLabel
(
    const word& description,
    const label computed,
    const label expected
)
{
    if (computed != expected)
    {
        FatalErrorInFunction
            << description << " does not match: " << endl
            << "  Expected value: " << expected << endl
            << "  Computed value: " << computed << endl
            << exit(FatalError);
    }
}

void checkBool
(
    const word& description,
    const bool computed,
    const bool expected
)
{
    if (computed != expected)
    {
        FatalErrorInFunction
            << description << " does not match: " << endl
            << "  Expected value: " << expected << endl
            << "  Computed value: " << computed << endl
            << exit(FatalError);
    }
}


// Tests univariate (1-D) access through the variadic and labelList
// overloads of operator(), both const and non-const, and confirms they
// agree with the underlying List storage order.
void testUnivariateAccess()
{
    Info<< "Testing univariate access" << endl;

    labelListList orders(4);
    orders[0] = labelList({0});
    orders[1] = labelList({1});
    orders[2] = labelList({2});
    orders[3] = labelList({3});

    mappedList<scalar> ml(4, orders, 0.0);

    ml(0) = 10.0;
    ml(1) = 20.0;
    ml(orders[2]) = 30.0;
    ml(labelList({3})) = 40.0;

    const mappedList<scalar>& cml = ml;

    checkScalar("ml(0)", ml(0), 10.0);
    checkScalar("ml(1)", ml(1), 20.0);
    checkScalar("ml(2)", ml(2), 30.0);
    checkScalar("ml(3)", ml(3), 40.0);

    checkScalar("cml(0)", cml(0), 10.0);
    checkScalar("cml(3)", cml(3), 40.0);

    checkScalar("ml(orders[0])", ml(orders[0]), 10.0);
    checkScalar("cml(orders[3])", cml(orders[3]), 40.0);

    // Underlying List storage order matches the constructor's index order
    checkScalar("ml[0]", ml[0], 10.0);
    checkScalar("ml[3]", ml[3], 40.0);

    Info<< "OK" << endl;
}


// Tests multivariate (2-D) access, verifying that the variadic and
// labelList overloads agree with each other. The chosen orders include
// (1 0), whose base-10 key (10) is what lets nDimensions be recovered
// correctly wherever it must be inferred from keys alone.
void testMultivariateAccess()
{
    Info<< "Testing multivariate access" << endl;

    labelListList orders(6);
    orders[0] = labelList({0, 0});
    orders[1] = labelList({1, 0});
    orders[2] = labelList({0, 1});
    orders[3] = labelList({1, 1});
    orders[4] = labelList({2, 0});
    orders[5] = labelList({0, 2});

    mappedList<scalar> ml(6, orders, 0.0);

    forAll(orders, oi)
    {
        ml(orders[oi]) = scalar(oi + 1);
    }

    checkScalar("ml(0,0)", ml(0, 0), 1.0);
    checkScalar("ml(1,0)", ml(1, 0), 2.0);
    checkScalar("ml(0,1)", ml(0, 1), 3.0);
    checkScalar("ml(1,1)", ml(1, 1), 4.0);
    checkScalar("ml(2,0)", ml(2, 0), 5.0);
    checkScalar("ml(0,2)", ml(0, 2), 6.0);

    const mappedList<scalar>& cml = ml;

    checkScalar("cml(1,0)", cml(1, 0), 2.0);
    checkScalar("cml(0,1)", cml(0, 1), 3.0);

    checkScalar
    (
        "ml(labelList{1,0}) vs ml(1,0)",
        ml(labelList({1, 0})),
        ml(1, 0)
    );

    Info<< "OK" << endl;
}


// Tests found(), including entries that are present, entries that are
// simply absent, and lists longer than the set's dimensionality (which
// must short-circuit to false rather than being looked up).
void testFound()
{
    Info<< "Testing found()" << endl;

    labelListList orders(3);
    orders[0] = labelList({0, 0});
    orders[1] = labelList({1, 0});
    orders[2] = labelList({0, 1});

    mappedList<scalar> ml(3, orders, 0.0);

    checkBool("found(0,0)", ml.found(0, 0), true);
    checkBool("found(1,0)", ml.found(1, 0), true);
    checkBool("found(0,1)", ml.found(0, 1), true);
    checkBool("found(1,1)", ml.found(1, 1), false);
    checkBool("found({1,0})", ml.found(labelList({1, 0})), true);
    checkBool("found({0,0,0})", ml.found(labelList({0, 0, 0})), false);

    Info<< "OK" << endl;
}


// Tests that setSize()/resize() rebuild the map while preserving the
// values already held at overlapping list positions.
void testSetSizeResize()
{
    Info<< "Testing setSize()/resize()" << endl;

    labelListList orders(2);
    orders[0] = labelList({0});
    orders[1] = labelList({1});

    mappedList<scalar> ml(2, orders, 0.0);
    ml(0) = 100.0;
    ml(1) = 200.0;

    labelListList newOrders(3);
    newOrders[0] = labelList({0});
    newOrders[1] = labelList({1});
    newOrders[2] = labelList({2});

    ml.setSize(3, newOrders);
    ml(2) = 300.0;

    checkScalar("ml(0) after setSize", ml(0), 100.0);
    checkScalar("ml(1) after setSize", ml(1), 200.0);
    checkScalar("ml(2) after setSize", ml(2), 300.0);
    checkBool("found(2) after setSize", ml.found(2), true);
    checkBool("found(5) after setSize", ml.found(5), false);

    labelListList resizeOrders(4);
    resizeOrders[0] = labelList({0});
    resizeOrders[1] = labelList({1});
    resizeOrders[2] = labelList({2});
    resizeOrders[3] = labelList({3});

    ml.resize(4, resizeOrders);
    ml(3) = 400.0;

    checkScalar("ml(3) after resize", ml(3), 400.0);

    Info<< "OK" << endl;
}


// Tests the indexes() accessor and the static listToWord() helper used to
// build the field-name suffix of an order list.
void testIndexesAndStaticHelpers()
{
    Info<< "Testing indexes() and static helpers" << endl;

    labelListList orders(3);
    orders[0] = labelList({0, 0});
    orders[1] = labelList({1, 0});
    orders[2] = labelList({0, 1});

    mappedList<scalar> ml(3, orders, 0.0);

    checkLabel("indexes().size()", ml.indexes().size(), 3);
    checkLabel("nDimensions()", ml.nDimensions(), 2);
    checkBool("indexes()[2] == (0 1)", ml.indexes()[2] == orders[2], true);

    if (mappedList<scalar>::listToWord(labelList({1, 2, 3})) != "123")
    {
        FatalErrorInFunction
            << "listToWord({1, 2, 3}) does not match the expected word "
            << "\"123\"" << endl
            << exit(FatalError);
    }

    Info<< "OK" << endl;
}


// Orders none of which has a non-zero leading component. The keys used to
// be packed a decimal digit to a dimension, and a count of dimensions read
// back from them came out short; the orders are now the keys.
void testLeadingZeroOrders()
{
    Info<< "Testing orders with leading zeros" << endl;

    labelListList orders(4);
    orders[0] = labelList({0, 0, 0});
    orders[1] = labelList({0, 1, 0});
    orders[2] = labelList({0, 1, 1});
    orders[3] = labelList({0, 2, 0});

    mappedList<scalar> ml(4, orders, 0.0);

    forAll(orders, oi)
    {
        ml(orders[oi]) = scalar(oi + 1);
    }

    checkScalar("ml(0,0,0)", ml(0, 0, 0), 1.0);
    checkScalar("ml(0,1,0)", ml(0, 1, 0), 2.0);
    checkScalar("ml(0,1,1)", ml(0, 1, 1), 3.0);
    checkScalar("ml(0,2,0)", ml(0, 2, 0), 4.0);

    Info<< "OK" << endl;
}


// Orders of ten or more, which the decimal packing of the keys could not
// tell apart: (1 10) and (11 0) were both key 110, and (0 10) was the key of
// (1 0). They are distinct elements now, found by their orders exactly.
void testOrdersOfTenOrMore()
{
    Info<< "Testing orders of ten or more" << endl;

    labelListList orders(5);
    orders[0] = labelList({1, 10});
    orders[1] = labelList({11, 0});
    orders[2] = labelList({1, 0});
    orders[3] = labelList({0, 10});
    orders[4] = labelList({255, 300});

    mappedList<scalar> ml(5, orders, 0.0);

    forAll(orders, oi)
    {
        ml(orders[oi]) = scalar(oi + 1);
    }

    checkScalar("ml(1,10)", ml(1, 10), 1.0);
    checkScalar("ml(11,0)", ml(11, 0), 2.0);
    checkScalar("ml(1,0)", ml(1, 0), 3.0);
    checkScalar("ml(0,10)", ml(0, 10), 4.0);
    checkScalar("ml(255,300)", ml(255, 300), 5.0);
    checkBool("found(11)", ml.found(11), true);
    checkBool("found(10,1)", ml.found(10, 1), false);
    checkBool("found(0,1)", ml.found(0, 1), false);
    checkBool
    (
        "found(1024,0), an order past what the key holds",
        ml.found(1024, 0),
        false
    );
    checkBool("found(-1,0)", ml.found(-1, 0), false);

    Info<< "OK" << endl;
}


// A list of more dimensions than a label has decimal digits, as the
// moments of a mixture fraction and a dozen reacting scalars are: the
// decimal keys overflowed there. With the mixture fraction up to order 40
// the orders no longer fit into a key of 64 bits either, and the list finds
// its elements by a hash of the orders and the orders themselves.
void testManyDimensions()
{
    Info<< "Testing a list of thirteen dimensions" << endl;

    const label nDims = 13;

    labelListList orders(nDims + 2);
    orders[0] = labelList(nDims, Zero);

    for (label dimi = 0; dimi < nDims; dimi++)
    {
        orders[dimi + 1] = labelList(nDims, Zero);
        orders[dimi + 1][dimi] = 1;
    }

    orders[nDims + 1] = labelList(nDims, Zero);
    orders[nDims + 1][0] = 40;

    mappedList<scalar> ml(nDims + 2, orders, 0.0);

    forAll(orders, oi)
    {
        ml(orders[oi]) = scalar(oi);
    }

    labelList last(nDims, Zero);
    last[nDims - 1] = 1;

    checkScalar("the last first-order moment", ml(last), scalar(nDims));
    checkScalar("the first one, by its leading order", ml(0, 1), 2.0);
    checkScalar("the zero-order moment, as ml(0)", ml(0), 0.0);

    labelList two(nDims, Zero);
    two[nDims - 1] = 2;

    checkBool("a second order that is not there", ml.found(two), false);
    checkScalar
    (
        "the order 40 of the first dimension",
        ml(40),
        scalar(nDims + 1)
    );
    checkBool("an order 41 that is not there", ml.found(41), false);

    Info<< "OK" << endl;
}


int main(int argc, char *argv[])
{
    Info<< "Testing mappedList\n" << endl;

    testUnivariateAccess();
    testMultivariateAccess();
    testFound();
    testSetSizeResize();
    testIndexesAndStaticHelpers();
    testLeadingZeroOrders();
    testOrdersOfTenOrMore();
    testManyDimensions();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
