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

Description
    Assertions shared by the tests of the multivariate moment inversions.

    The check that carries the weight is moment conservation: a multivariate
    inversion builds a quadrature from the moments of its own moment set, so
    it has to reproduce every one of them. That is a property of the method
    rather than a value recorded from a previous run.

    It only holds for the moments the method controls. The moment sets the
    tests build from a quadrature of their own carry moments outside that
    set, and the value an inversion returns for one of those is the closure
    it applies, not an error: CHyQMOM returns moment (1 0 2) of a random
    twenty-seven node quadrature with a relative error of two hundred, and
    is behaving as designed. Only the moments of the method are asserted.

\*---------------------------------------------------------------------------*/

#ifndef multivariateMomentTest_H
#define multivariateMomentTest_H

#include "labelList.H"
#include "scalarList.H"
#include "messageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Name of a moment of the given order, as moment.012
inline word momentName(const labelList& momentOrder)
{
    word name("moment.");

    forAll(momentOrder, dimi)
    {
        name += Foam::name(momentOrder[dimi]);
    }

    return name;
}


//- Report a computed value against the expected one, and stop with a fatal
//  error if they differ by more than the tolerance, relative to the
//  magnitude of the expected value
inline void compareScalar
(
    const scalar computed,
    const scalar expected,
    const scalar tolerance,
    const string& name
)
{
    const scalar difference = mag(computed - expected);
    const scalar bound = tolerance*max(mag(expected), SMALL);

    Info<< "  " << name << " = " << computed
        << ", expected " << expected << endl;

    if (difference > bound)
    {
        FatalErrorInFunction
            << "The value of " << name << " is not the expected one." << nl
            << "    Computed: " << computed << nl
            << "    Expected: " << expected << nl
            << "    Difference: " << difference << nl
            << "    Tolerance: " << bound << nl
            << exit(FatalError);
    }
}


//- Assert that the quadrature reproduces every moment of momentOrders.
//  reconstructed and expected are read with the order as the key, so the
//  order the moments are listed in does not matter.
//
//  The error is measured against the largest moment of the set rather than
//  against each moment on its own, because a moment set is free to contain
//  moments that are exactly zero, and an uncorrelated distribution does:
//  a purely relative tolerance would ask those to be reproduced to no error
//  at all. The criterion is that every moment is reproduced to within the
//  tolerance times the magnitude of the largest moment of the set.
//  knownFailures lists the moments a method is known not to conserve. They
//  are reported and not asserted, so that a defect waiting to be fixed does
//  not hide every other moment behind it. One that starts conserving is
//  reported too, because the pin has to be removed once it is fixed: left
//  in place it would hide the defect coming back.
template<class reconstructedType, class expectedType>
void checkMomentConservation
(
    const reconstructedType& reconstructed,
    const expectedType& expected,
    const labelListList& momentOrders,
    const scalar tolerance,
    const word& what,
    const labelListList& knownFailures = labelListList()
)
{
    Info<< "\nVerifying moment conservation of " << what << endl;

    scalar scale = 0.0;

    forAll(momentOrders, mi)
    {
        scale = max(scale, mag(expected(momentOrders[mi])));
    }

    const scalar bound = tolerance*max(scale, SMALL);

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        const scalar computed = reconstructed(momentOrder);
        const scalar reference = expected(momentOrder);
        const scalar difference = mag(computed - reference);
        const bool conserved = difference <= bound;

        bool known = false;

        forAll(knownFailures, ki)
        {
            if (knownFailures[ki] == momentOrder)
            {
                known = true;
                break;
            }
        }

        Info<< "  " << momentName(momentOrder) << " = " << computed
            << ", expected " << reference
            << (known ? (conserved ? "  (pinned)" : "  (known failure)") : "")
            << endl;

        if (known)
        {
            if (conserved)
            {
                WarningInFunction
                    << momentName(momentOrder) << " of " << what
                    << " is pinned as a known failure, and is now"
                    << " conserved." << nl
                    << "    Remove it from the known failures of the test,"
                    << " so that the defect is caught if it returns." << nl
                    << endl;
            }

            continue;
        }

        if (!conserved)
        {
            FatalErrorInFunction
                << "The inversion does not conserve "
                << momentName(momentOrder) << " of " << what << "." << nl
                << "    Computed: " << computed << nl
                << "    Expected: " << reference << nl
                << "    Difference: " << difference << nl
                << "    Tolerance: " << bound
                << " (" << tolerance << " of the largest moment, "
                << scale << ")" << nl
                << exit(FatalError);
        }
    }
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
