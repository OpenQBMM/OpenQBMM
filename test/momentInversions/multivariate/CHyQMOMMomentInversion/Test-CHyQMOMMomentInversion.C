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
    Test-CHyQMOMMomentInversion

Description
    Test the conditional hyperbolic moment inversion classes.

    A quadrature of twenty-seven nodes is drawn, its moments are computed,
    and both CHyQMOM and CHyQMOM+ are asked to invert them. Each has to
    reproduce the moments of its own moment set.

    The quadrature carries more moments than either method controls, and
    the value an inversion returns for one of the others is the closure it
    applies rather than an error, so those are reported and not asserted.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "mappedLists.H"
#include "CHyQMOMMomentInversion.H"
#include "CHyQMOMPlusMomentInversion.H"
#include "Random.H"
#include <cfenv>
#include "multivariateMomentTest.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Invert a moment set of very large moments and check that it raises no
//  floating point exception.
//
//  The set is the one the size conditioned inversion of the bidisperse
//  granular cooling case hands over, of zero-order moment one and second
//  and fourth order moments of 1e102 and 1e206. Testing realizability as
//  c2 c4 < c2^3 + c3^2 forms a product that is not a finite double, and the
//  run of that case ended on the exception. The condition is the same one
//  written on the standardised moments, eta < 1 + q^2, and nothing there
//  leaves the range of a double.
//
//  The check is on the exception flags rather than on the moments that come
//  back: an inversion that overflows still returns numbers, and it is the
//  overflow itself that stops a solver, since OpenFOAM traps it.
void testExtremeMomentSet(const dictionary& dict)
{
    Info<< "\n\nInverting a moment set of very large moments" << endl;

    const labelListList momentOrders
    (
        multivariateMomentInversions::CHyQMOM::getMomentOrders(3)
    );

    const labelListList nodeIndexes
    (
        multivariateMomentInversions::CHyQMOM::getNodeIndexes(3)
    );

    const scalarList values
    ({
        1.0,
       -1.3784e+35,   1.08747e+35, -1.52797e+36,
        3.63077e+102, 4.00183e+100, 4.00183e+100,
        3.63077e+102, 4.00183e+100, 3.63077e+102,
       -1.15483e+139, 1.87662e+139, -1.67508e+139,
        5.04075e+205, 5.04075e+205, 5.04075e+205
    });

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
        moments(momentOrders[mi]) = values[mi];
    }

    multivariateMomentInversions::CHyQMOM inverter
    (
        dict, momentOrders, nodeIndexes, {0, 1, 2}
    );

    std::feclearexcept(FE_ALL_EXCEPT);

    inverter.invert(moments);

    const int raised =
        std::fetestexcept(FE_OVERFLOW | FE_DIVBYZERO | FE_INVALID);

    if (raised)
    {
        FatalErrorInFunction
            << "The inversion of a moment set of very large moments raised "
            << "a floating point exception." << nl
            << "    Overflow: " << bool(raised & FE_OVERFLOW) << nl
            << "    Division by zero: " << bool(raised & FE_DIVBYZERO) << nl
            << "    Invalid: " << bool(raised & FE_INVALID) << nl
            << exit(FatalError);
    }

    Info<< "  no floating point exception was raised" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Invert the moments with the given method and assert that it conserves
//  every moment of its own moment set
template<class inversionType>
void testInversion
(
    const word& what,
    const dictionary& dict,
    const multivariateMomentSet& moments,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes,
    const label nDims,
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

            for (label dimi = 0; dimi < nDims; dimi++)
            {
                cmpt *= pow(abscissae(nodeIndex)[dimi], momentOrder[dimi]);
            }

            m += cmpt;
        }

        reconstructed(momentOrder) = m;
    }

    // The moments of the method, which it has to conserve
    const labelListList controlled(inversionType::getMomentOrders(nDims));

    checkMomentConservation
    (
        reconstructed,
        moments,
        controlled,
        tolerance,
        what,
        knownFailures
    );

    // The remaining moments of the quadrature are closed by the method
    Info<< "\nMoments outside the set of " << what
        << ", closed rather than conserved:" << endl;

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        bool isControlled = false;

        forAll(controlled, ci)
        {
            if (controlled[ci] == momentOrder)
            {
                isControlled = true;
                break;
            }
        }

        if (!isControlled)
        {
            Info<< "  " << momentName(momentOrder) << " = "
                << reconstructed(momentOrder)
                << ", of a quadrature with " << moments(momentOrder) << endl;
        }
    }
}


int main()
{
    #include "createFields.H"

    mappedList<scalarList> x(nNodes, nodeIndexes, scalarField(nDims, Zero));
    mappedList<scalar> w(nNodes, nodeIndexes, 0.0);

    // A fixed seed, so that the moments the inversions are asked for are the
    // same on every run and on every machine
    Random rndGen(20260904);

    forAll(x, nodei)
    {
        w[nodei] = rndGen.sample01<scalar>();

        forAll(x[nodei], dimi)
        {
            x[nodei][dimi] = 2.0*rndGen.sample01<scalar>() - 1.0;
        }
    }

    multivariateMomentSet moments
    (
        nMoments,
        momentOrders,
        List<supportType>(momentOrders[0].size(), supportType::R),
        SMALL,
        SMALL
    );

    Info<< "Moments of the quadrature:" << endl;

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scalar m = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            scalar cmpt = w(nodeIndex);

            forAll(nodeIndex, dimi)
            {
                cmpt *= pow(x(nodeIndex)[dimi], momentOrder[dimi]);
            }

            m += cmpt;
        }

        moments(momentOrder) = m;

        Info<< "  " << momentName(momentOrder) << ": " << m << endl;
    }

    const scalar tolerance = 1e-10;

    testInversion<multivariateMomentInversions::CHyQMOM>
    (
        "CHyQMOM",
        quadratureProperties,
        moments,
        momentOrders,
        nodeIndexes,
        velocityIndexes,
        nDims,
        tolerance
    );

    testInversion<multivariateMomentInversions::CHyQMOMPlus>
    (
        "CHyQMOMPlus",
        quadratureProperties,
        moments,
        momentOrders,
        nodeIndexes,
        velocityIndexes,
        nDims,
        tolerance
    );

    testExtremeMomentSet(quadratureProperties);

    Info<< "\n\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
