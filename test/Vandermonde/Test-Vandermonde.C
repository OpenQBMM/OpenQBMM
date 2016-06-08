/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 Alberto Passalacqua
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
    Test-ExtendedMomentInversion.C

Description
    Test the extendedMomentInversion class and its subclasses.

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
    label n = 5;
    scalarDiagonalMatrix A(n);
    forAll(A, i)
    {
        A[i] = scalar(rand())/scalar(RAND_MAX);
    }

    Vandermonde vm(A);
    Vandermonde V(vm());

    Info<< nl << "Initial Vector: " << A << endl;
    Info<< "Vector constructed from square Vandermonde system: " << nl << "\t"
        << V << endl;

    scalarSquareMatrix invVM = vm.inv();

    Info<< nl
        << "Vandermonde system: " << endl;

    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            Info<< vm(i,j) << ", \t";
        }
        Info<< endl;
    }

    Info<< nl
        << "Inverted Vandermonde matrix:" << endl;

    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            Info<< invVM(i,j) << ", \t";
        }
        Info<< endl;
    }

    scalarRectangularMatrix VM = SVDinv(invVM);

    scalar error = 0.0;
    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            error += sqr(VM(i,j) - vm(i,j));
        }
    }

    error = Foam::sqrt(error);
    Info<< nl << "Total magnitude of error when Vandermonde system is " << nl
        << "inverted twice: " << error << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
