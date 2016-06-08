/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/

#include "Vandermonde.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Vandermonde::Vandermonde
(
    const scalarDiagonalMatrix& A
)
:
    scalarDiagonalMatrix(A)
{}


Foam::Vandermonde::Vandermonde
(
    const scalarSquareMatrix& A
)
:
    scalarDiagonalMatrix(A.m())
{
    if (A.m() != A.n())
    {
        FatalErrorInFunction
            << "Refrence matrix is not a Vandermonde system." << nl
            << "Matrix is not square."
            << abort(FatalError);
    }

    scalar sumRow = 0.0;
    scalar sumCol = 0.0;
    for (label i = 0; i < size(); i++)
    {
        sumRow += A(0,i);
        sumCol += A(i,0);
    }

    if (sumRow == A.m())
    {
        for (label i = 0; i < size(); i++)
        {
            (*this)[i] = A(1,i);
        }
    }
    else if (sumCol == A.m())
    {
        for (label i = 0; i < size(); i++)
        {
            (*this)[i] = A(i,1);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Refrence matrix is not a Vandermonde system." << nl
            << "First row is not composed of ones."
            << abort(FatalError);
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Vandermonde::~Vandermonde()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::Vandermonde::solve
(
    scalarDiagonalMatrix& x,
    const scalarDiagonalMatrix& source
)
{
    const label n = size();
    scalarDiagonalMatrix c(n, 0.0);

    scalar t = 1.0;
    scalar r = 1.0;
    scalar s = 0.0;
    scalar xx = 0.0;

    if (n == 1)
    {
        x[0] = source[0];
    }
    else
    {
        c[n-1] = -(*this)[0];

        for (label i = 1; i < n; i++)
        {
            xx = -(*this)[i];

            for (label j = (n-i-1); j < (n-1); j++)
            {
                c[j] += xx*c[j+1];
            }
            c[n-1] += xx;
        }

        for (label i = 0; i < n; i++)
        {
            xx = (*this)[i];
            t = 1.0;
            r = 1.0;
            s = source[n-1];

            for (label j = n - 1; j > 0; j--)
            {
                r = c[j] + r*xx;
                s += r*source[j-1];
                t = r + t*xx;
            }

            x[i] = s/t;
        }
    }
}

Foam::scalarSquareMatrix Foam::Vandermonde::inv()
{
    const label n = size();

    scalarSquareMatrix inv(n);
    scalarDiagonalMatrix source(n);
    scalarDiagonalMatrix x(n);

    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            if (i != j)
            {
                source[j] = 0.0;
            }
            else
            {
                source[j] = 1.0;
            }
        }

        solve(x, source);

        for (label j = 0; j < n; j++)
        {
            inv[j][i] = x[j];
        }
    }

    return inv;
}


// ************************************************************************* //
