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
    Copyright (C) 2019-2023 Alberto Passalacqua
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
    const scalar expectedError = 8.53672259787214e-13;
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
    Vandermonde V(Vm());

    Info<< "Initial vector: " << A << endl;

    Info<< "Vector constructed from square Vandermonde matrix: " << V << endl;

    scalarSquareMatrix invVm = Vm.inv();

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

    if (mag(error - expectedError) >= SMALL)
    {
        FatalErrorInFunction
                << "The error accumulated during two inversions is large: " 
                << endl
                << "  Expected maximum error: " << expectedError
                << endl
                << "  Actual error: " << error
                << endl
                << exit(FatalError);
    }

    Info<< nl << "Total magnitude of error when Vandermonde matrix is"
        << nl << "inverted twice: " << error << ", OK." << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
