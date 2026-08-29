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
    Test-MappedPtrList

Description
    Test the mappedPtrList class.

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "mappedPtrList.H"

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


// Tests construction from (size, indexes), population through the
// raw-pointer set() overload, and access through the variadic and
// labelList overloads of operator(), both const and non-const.
void testConstructAndAccess()
{
    Info<< "Testing construction and access" << endl;

    labelListList orders(6);
    orders[0] = labelList({0, 0});
    orders[1] = labelList({1, 0});
    orders[2] = labelList({0, 1});
    orders[3] = labelList({1, 1});
    orders[4] = labelList({2, 0});
    orders[5] = labelList({0, 2});

    mappedPtrList<scalar> mp(6, orders);

    forAll(orders, oi)
    {
        mp.set(orders[oi], new scalar(oi + 1));
    }

    checkScalar("mp(0,0)", mp(0, 0), 1.0);
    checkScalar("mp(1,0)", mp(1, 0), 2.0);
    checkScalar("mp(0,1)", mp(0, 1), 3.0);
    checkScalar("mp(1,1)", mp(1, 1), 4.0);
    checkScalar("mp(2,0)", mp(2, 0), 5.0);
    checkScalar("mp(0,2)", mp(0, 2), 6.0);

    const mappedPtrList<scalar>& cmp = mp;
    checkScalar("cmp(1,0)", cmp(1, 0), 2.0);

    checkScalar
    (
        "mp(labelList{2,0}) vs mp(2,0)",
        mp(labelList({2, 0})),
        mp(2, 0)
    );

    checkBool("set(0)", mp.set(label(0)), true);
    checkBool("set({0,0})", mp.set(orders[0]), true);

    Info<< "OK" << endl;
}


// Tests the autoPtr overload of set(). This previously bound OpenFOAM's
// deprecated set(label, autoPtr<T>&) overload instead of the intended
// set(label, autoPtr<T>&&) one, because the named autoPtr parameter was
// passed on as an lvalue.
void testSetAutoPtr()
{
    Info<< "Testing set() with autoPtr" << endl;

    labelListList orders(2);
    orders[0] = labelList({0});
    orders[1] = labelList({1});

    mappedPtrList<scalar> mp(2, orders);

    autoPtr<scalar> value(new scalar(42.0));
    mp.set(orders[1], std::move(value));

    checkBool("set(1) after autoPtr set", mp.set(label(1)), true);
    checkScalar("mp(1) after autoPtr set", mp(1), 42.0);

    Info<< "OK" << endl;
}


// Tests found(), analogous to mappedList's.
void testFound()
{
    Info<< "Testing found()" << endl;

    labelListList orders(3);
    orders[0] = labelList({0, 0});
    orders[1] = labelList({1, 0});
    orders[2] = labelList({0, 1});

    mappedPtrList<scalar> mp(3, orders);

    checkBool("found(0,0)", mp.found(0, 0), true);
    checkBool("found(1,0)", mp.found(1, 0), true);
    checkBool("found(1,1)", mp.found(1, 1), false);
    checkBool("found({0,0,0})", mp.found(labelList({0, 0, 0})), false);

    Info<< "OK" << endl;
}


// Tests setMap(const labelListList&), the path quadratureApproximation
// uses to key its node list after building it from a dictionary Istream:
// re-keys an existing list with a differently-shaped set of orders and
// confirms the new map, not the old one, is what found() sees.
void testSetMapFromIndexes()
{
    Info<< "Testing setMap(const labelListList&)" << endl;

    labelListList orders(3);
    orders[0] = labelList({0});
    orders[1] = labelList({1});
    orders[2] = labelList({2});

    mappedPtrList<scalar> mp(3, orders);

    labelListList newOrders(3);
    newOrders[0] = labelList({0, 0});
    newOrders[1] = labelList({1, 0});
    newOrders[2] = labelList({0, 1});

    mp.setMap(newOrders);

    checkBool("found(0,0) after setMap(indexes)", mp.found(0, 0), true);
    checkBool("found(1,0) after setMap(indexes)", mp.found(1, 0), true);
    checkBool("found(0,1) after setMap(indexes)", mp.found(0, 1), true);
    checkBool("found(2,2) after setMap(indexes)", mp.found(2, 2), false);

    Info<< "OK" << endl;
}


// Tests setMap(const Map<label>&), the path used when a new mappedPtrList
// reuses the key map already computed for another one (e.g. moment
// advection reusing a quadrature's node map). Also confirms that a
// second call with a lower-dimensional map does not retain the
// dimensionality left over from the first call: without resetting
// nDimensions before recomputing it, a query for order {2} would be
// packed using the stale dimensionality and incorrectly miss the map.
void testSetMapFromMap()
{
    Info<< "Testing setMap(const Map<label>&)" << endl;

    labelListList orders2D(3);
    orders2D[0] = labelList({0, 0});
    orders2D[1] = labelList({1, 0});
    orders2D[2] = labelList({0, 1});

    mappedPtrList<scalar> source2D(3, orders2D);

    mappedPtrList<scalar> mp(3, orders2D);
    mp.setMap(source2D.map());

    checkBool("found(1,0) after setMap(2-D map)", mp.found(1, 0), true);
    checkBool("found(0,1) after setMap(2-D map)", mp.found(0, 1), true);

    labelListList orders1D(3);
    orders1D[0] = labelList({0});
    orders1D[1] = labelList({1});
    orders1D[2] = labelList({2});

    mappedPtrList<scalar> source1D(3, orders1D);
    mp.setMap(source1D.map());

    checkBool
    (
        "found({2}) after re-setMap to a 1-D map",
        mp.found(labelList({2})),
        true
    );

    Info<< "OK" << endl;
}


// Tests setSize()/resize().
void testSetSizeResize()
{
    Info<< "Testing setSize()/resize()" << endl;

    labelListList orders(2);
    orders[0] = labelList({0});
    orders[1] = labelList({1});

    mappedPtrList<scalar> mp(2, orders);
    mp.set(orders[0], new scalar(100.0));
    mp.set(orders[1], new scalar(200.0));

    labelListList newOrders(3);
    newOrders[0] = labelList({0});
    newOrders[1] = labelList({1});
    newOrders[2] = labelList({2});

    mp.setSize(3, newOrders);
    mp.set(newOrders[2], new scalar(300.0));

    checkScalar("mp(0) after setSize", mp(0), 100.0);
    checkScalar("mp(1) after setSize", mp(1), 200.0);
    checkScalar("mp(2) after setSize", mp(2), 300.0);

    Info<< "OK" << endl;
}


int main(int argc, char *argv[])
{
    Info<< "Testing mappedPtrList\n" << endl;

    testConstructAndAccess();
    testSetAutoPtr();
    testFound();
    testSetMapFromIndexes();
    testSetMapFromMap();
    testSetSizeResize();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
