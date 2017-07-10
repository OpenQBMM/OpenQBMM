/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 Alberto Passalacqua
     \\/     M anipulation  |
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
    Test-EigenSolver

Description
    Example application of the eigenSolver class.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "scalarMatrices.H"
#include "eigenSolver.H"

int main(int argc, char *argv[])
{
    scalarSquareMatrix A(3, scalar(0));

    A[0][0] = 1.0;
    A[0][1] = 2.0;
    A[0][2] = 3.0;

    A[1][0] = 2.0;
    A[1][1] = 7.0;
    A[1][2] = 4.0;

    A[2][0] = 3.0;
    A[2][1] = 4.0;
    A[2][2] = 9.0;

    Info<< "Input matrix" << endl << A << endl;

    // Since the matrix is not symmtric, there are two possible ways of
    // constructing the eigenSolver:
    //
    // a) Using the general constructor, with symmetry check

    eigenSolver eig(A);

    // b) Skypping the symmetry check by specifying the matrix is not
    //    symmetric
    //
    // eigenSolver eig(A, false);


    // The imaginary part of eigenvalues can be retrieved with
    //
    // scalarDiagonalMatrix eigenIm(eig.eigenvaluesIm());
    //
    // while eigenvectors are retrived with
    //
    // scalarSquareMatrix eigenvec(eig.eigenvectors());

    // Printing the real part of eigenvalues
    Info<< "\nReal part of eigenvalues: " << endl << eig.eigenvaluesRe()
        << endl;

    Info<< "\nImaginary part of eigenvalues: " << endl
        << eig.eigenvaluesIm() << endl;

    Info<< "\nEigenvectors" << endl << eig.eigenvectors() << "\n\n";

    Info<< "End\n" << endl;

    return 0;
}
