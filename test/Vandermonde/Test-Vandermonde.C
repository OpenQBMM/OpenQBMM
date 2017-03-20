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
    label n = 5;
    scalarDiagonalMatrix A(n);

    forAll(A, i)
    {
        A[i] = scalar(rand())/scalar(RAND_MAX);
    }

    Vandermonde Vm(A);
    Vandermonde V(Vm());

    Info<< nl << "Initial vector: " << A << endl;

    Info<< "Vector constructed from square Vandermonde matrix: " << nl << "\t"
        << V << endl;

    scalarSquareMatrix invVm = Vm.inv();

    Info<< nl << "Vandermonde matrix: " << endl;

    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            Info<< Vm(i, j) << ", \t";
        }

        Info<< endl;
    }

    Info<< nl << "Inverted Vandermonde matrix:" << endl;

    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            Info<< invVm(i, j) << ", \t";
        }

        Info<< endl;
    }

    scalarRectangularMatrix svdInv = SVDinv(invVm);

    scalar error = 0.0;

    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            error += sqr(svdInv(i, j) - Vm(i, j));
        }
    }

    error = Foam::sqrt(error);

    Info<< nl << "Total magnitude of error when Vandermonde system is "
        << nl << "inverted twice: " << error << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
