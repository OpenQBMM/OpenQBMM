/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2026 Alberto Passalacqua
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
    Test-Vandermonde.C

Description
    Test the Vandermonde class.

\*---------------------------------------------------------------------------*/

#include <cmath>
#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "Random.H"
#include "Vandermonde.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<< setprecision(16);

    Info<< "Testing Vandermonde\n" << endl;
    Info<< "-------------------\n" << endl;

    label matrixSize = 5;

    scalarSquareMatrix expectedInverse(matrixSize);

    expectedInverse(0, 0) = 1;
    expectedInverse(0, 1) = -2.083333333333333;
    expectedInverse(0, 2) = 1.458333333333333;
    expectedInverse(0, 3) = -0.4166666666666667;
    expectedInverse(0, 4) = 0.04166666666666666;

    expectedInverse(1, 0) = 0.0;
    expectedInverse(1, 1) = 4.0;
    expectedInverse(1, 2) = -4.333333333333333;
    expectedInverse(1, 3) = 1.5;
    expectedInverse(1, 4) = -0.1666666666666667;

    expectedInverse(2, 0) = 0.0;
    expectedInverse(2, 1) = -3.0;
    expectedInverse(2, 2) = 4.75;
    expectedInverse(2, 3) = -2.0;
    expectedInverse(2, 4) = 0.25;

    expectedInverse(3, 0) = 0.0;
    expectedInverse(3, 1) = 1.333333333333333;
    expectedInverse(3, 2) = -2.333333333333333;
    expectedInverse(3, 3) = 1.166666666666667;
    expectedInverse(3, 4) = -0.1666666666666667;

    expectedInverse(4, 0) = 0.0;
    expectedInverse(4, 1) = -0.25;
    expectedInverse(4, 2) = 0.4583333333333333;
    expectedInverse(4, 3) = -0.25;
    expectedInverse(4, 4) = 0.04166666666666666;

    scalarDiagonalMatrix A(matrixSize);

    forAll(A, i)
    {
        A[i] = i;
    }

    Vandermonde Vm(A);
    Vandermonde V(Vm(), true);

    // solve() and invert() are const: exercise them through a const
    // reference so a regression on that is caught at compile time.
    const Vandermonde& constVm = Vm;

    Info<< "Initial vector: " << A << endl;

    Info<< "Vector constructed from square Vandermonde matrix: " << V << endl;

    scalarSquareMatrix invVm = constVm.invert();

    Info<< nl << "Vandermonde matrix:\n" << endl;

    for (label i = 0; i < matrixSize; i++)
    {
        for (label j = 0; j < matrixSize; j++)
        {
            Info<< "  " << Vm(i, j) << ", \t";
        }

        Info<< endl;
    }

    Info<< nl << "Inverse of Vandermonde matrix\n" << endl;

    for (label i = 0; i < matrixSize; i++)
    {
        for (label j = 0; j < matrixSize; j++)
        {
            Info<< "  " << invVm(i, j) << ", \t";
        }

        Info<< endl;
    }

    for(int i = 0; i < matrixSize; i++)
    {
        for(int j = 0; j < matrixSize; j++)
        {
            scalar magDiff = mag(invVm(i, j) - expectedInverse(i,j));

            if (magDiff > SMALL)
            {
                FatalErrorInFunction
                << "Element (" << i << ", " << j << ") in inverse matrix "
                << "differs from expected value: "
                << endl
                << "  Expected value: " << expectedInverse(i,j)
                << endl
                << "  Computed value: " << invVm(i, j)
                << endl
                << exit(FatalError);
            }
        }
    }

    Info<< endl << "Inverse matrix matches." << endl;

    scalarRectangularMatrix svdInv = SVDinv(invVm);

    scalar error = 0.0;

    for (label i = 0; i < matrixSize; i++)
    {
        for (label j = 0; j < matrixSize; j++)
        {
            error += sqr(svdInv(i, j) - Vm(i, j));
        }
    }

    error = Foam::sqrt(error);

    // The exact value of this round-trip error depends on the platform's
    // BLAS/LAPACK implementation, so check that it stays comfortably
    // small rather than pinning it to one machine's result to the last
    // digit.
    const scalar maxExpectedError = 1.0e-9;

    if (error > maxExpectedError)
    {
        FatalErrorInFunction
                << "The error accumulated during two inversions is too "
                << "large: "
                << endl
                << "  Maximum expected error: " << maxExpectedError
                << endl
                << "  Actual error: " << error
                << endl
                << exit(FatalError);
    }

    Info<< nl << "Total magnitude of error when Vandermonde matrix is"
        << nl << "inverted twice: " << error << ", OK." << endl;


    // Test the Vandermonde solve method
    scalarDiagonalMatrix source(matrixSize);
    scalarDiagonalMatrix x(matrixSize);

    for (label i = 0; i < matrixSize; i++)
    {
        source[i] = i;
    }

    constVm.solve(x, source);

    Info<< "\nTesting solve method: " << endl;

    // Print the known vector
    Info<< "\nKnown vector: " << source << endl;

    // Print the solution vector
    Info<< "\nSolution vector: " << x << endl;

    // Verify the solution by multiplying the matrix by the solution vector
    scalarDiagonalMatrix verification(matrixSize);
    
    for (label i = 0; i < matrixSize; i++)
    {
        verification[i] = 0.0;
        for (label j = 0; j < matrixSize; j++)
        {
            verification[i] += Vm(i, j) * x[j];
        }
    }

    Info<< "\nCalculated source vector: " << verification << endl;

    // Calculate the difference between the components of the source vector 
    /// and of the recomputed source vector (absolute value)
    scalarDiagonalMatrix difference(matrixSize);

    for (label i = 0; i < matrixSize; i++)
    {
        difference[i] = mag(verification[i] - source[i]);
    }

    // Print the difference vector
    Info<< "\nDifference between source vector components: "
        << difference << endl;

    // Test that a structurally valid Vandermonde matrix with large,
    // independently computed entries is accepted by isVandermonde(). An
    // absolute tolerance previously rejected this: matrix entries grow
    // like base^i, so round-off between std::pow() and the repeated
    // multiplication isVandermonde() checks against is not comparable to
    // a fixed tolerance once the abscissae move away from unit scale.
    Info<< "\nTesting isVandermonde() with large-magnitude entries: "
        << endl;

    const label largeSize = 6;

    scalarDiagonalMatrix largeBase(largeSize);
    largeBase[0] = 0.0;
    largeBase[1] = 1.0;
    largeBase[2] = 2.0;
    largeBase[3] = 5.0;
    largeBase[4] = 9.1;
    largeBase[5] = 12.4;

    scalarSquareMatrix largeA(largeSize);

    for (label j = 0; j < largeSize; j++)
    {
        for (label i = 0; i < largeSize; i++)
        {
            largeA(i, j) = std::pow(largeBase[j], scalar(i));
        }
    }

    // Constructing with checkVandermonde = true would FatalError if the
    // matrix were (incorrectly) rejected.
    Vandermonde largeV(largeA, true);

    Info<< "OK" << endl;

    // Test the degenerate 1x1 case: the constructor and isVandermonde()
    // previously read a non-existent second matrix row (A[1][...]) for a
    // 1x1 matrix, which is undefined behaviour.
    Info<< "\nTesting the degenerate 1x1 case: " << endl;

    scalarSquareMatrix A1(1);
    A1(0, 0) = 1.0;

    Vandermonde V1(A1, true);

    scalarDiagonalMatrix source1(1);
    scalarDiagonalMatrix x1(1);
    source1[0] = 7.0;

    V1.solve(x1, source1);

    if (mag(x1[0] - source1[0]) > SMALL)
    {
        FatalErrorInFunction
            << "Solving the 1x1 Vandermonde system did not return the "
            << "source value: "
            << endl
            << "  Expected value: " << source1[0]
            << endl
            << "  Computed value: " << x1[0]
            << endl
            << exit(FatalError);
    }

    Info<< "OK" << endl;

    // The conditional solve the quadratures share: from the moments of a
    // quantity weighted by the nodes of a quadrature, its value at each
    // node. Two of the four nodes carry no weight, and sit at the same
    // abscissa, which made the full system singular: they are left out,
    // and the other two are recovered from two moments.
    Info<< "\nTesting the conditional solve: " << endl;

    {
        const scalarList x({0.0, 0.5, 2.0, 0.0});
        const scalarList w({0.0, 3.0, 1.5, 0.0});
        const scalarList value({9.0, -1.25, 4.0, 9.0});

        // The moments of orders 0 and 1 of the quantity
        scalarList m(4, Zero);

        forAll(x, i)
        {
            m[0] += w[i]*value[i];
            m[1] += w[i]*x[i]*value[i];
        }

        labelList nodes;
        scalarSquareMatrix invVR;

        const bool solved =
            Vandermonde::conditionalInverse(x, w, 0, SMALL, nodes, invVR);

        if (!solved || nodes != labelList({1, 2}))
        {
            FatalErrorInFunction
                << "The conditional solve did not keep the two nodes that "
                << "carry weight: solved " << solved << ", nodes " << nodes
                << endl
                << exit(FatalError);
        }

        forAll(nodes, i)
        {
            scalar recovered = 0;

            forAll(nodes, k)
            {
                recovered += invVR(i, k)*m[k];
            }

            const scalar expected = value[nodes[i]];

            if (mag(recovered - expected) > 1.0e-12*mag(expected))
            {
                FatalErrorInFunction
                    << "The conditional solve recovered " << recovered
                    << " at node " << nodes[i] << ", where it is "
                    << expected << endl
                    << exit(FatalError);
            }
        }

        // A node that carries a little, below the weight its value is read
        // back at, still contributes to the moments. It stays in the system,
        // so the others are recovered exactly; left out, it would give them
        // its share, which its large abscissa makes a large one.
        {
            const scalarList xSmall({0.5, 2.0, 40.0});
            const scalarList wSmall({3.0, 1.5, 1.0e-10});
            const scalarList valueSmall({-1.25, 4.0, 7.0});

            scalarList mSmall(3, Zero);

            forAll(xSmall, i)
            {
                scalar xk = 1;

                for (label k = 0; k < 3; k++)
                {
                    mSmall[k] += wSmall[i]*xk*valueSmall[i];
                    xk *= xSmall[i];
                }
            }

            const bool solvedSmall =
                Vandermonde::conditionalInverse
                (
                    xSmall, wSmall, 0, 1.0e-6, nodes, invVR
                );

            if (!solvedSmall || nodes.size() != 3)
            {
                FatalErrorInFunction
                    << "The conditional solve left a node that carries a "
                    << "little out of the system: nodes " << nodes << endl
                    << exit(FatalError);
            }

            for (label i = 0; i < 2; i++)
            {
                scalar recovered = 0;

                for (label k = 0; k < 3; k++)
                {
                    recovered += invVR(i, k)*mSmall[k];
                }

                if (mag(recovered - valueSmall[i]) > 1.0e-10)
                {
                    FatalErrorInFunction
                        << "The conditional solve recovered " << recovered
                        << " at node " << i << ", where it is "
                        << valueSmall[i] << endl
                        << exit(FatalError);
                }
            }

            scalar lastRow = 0;

            for (label k = 0; k < 3; k++)
            {
                lastRow += mag(invVR(2, k));
            }

            if (lastRow != 0)
            {
                FatalErrorInFunction
                    << "The row of a node below the weight its value is "
                    << "read back at is not zero" << endl
                    << exit(FatalError);
            }
        }

        // Two nodes that carry weight at the same abscissa cannot be told
        // apart, and the solve says so rather than inverting a singular
        // system
        const scalarList xSame({0.5, 0.5, 2.0});
        const scalarList wSame({1.0, 2.0, 1.0});

        if
        (
            Vandermonde::conditionalInverse
            (
                xSame, wSame, 0, SMALL, nodes, invVR
            )
        )
        {
            FatalErrorInFunction
                << "The conditional solve accepted two nodes of weight at "
                << "the same abscissa" << endl
                << exit(FatalError);
        }
    }

    Info<< "OK" << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
