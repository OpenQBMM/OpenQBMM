/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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
    scalarDiagonalMatrix(A),
    n_(A.size())
{}


Foam::Vandermonde::Vandermonde
(
    const scalarSquareMatrix& A,
    const bool checkVandermonde
)
:
    scalarDiagonalMatrix(A.m()),
    n_(A.m())
{
    if (checkVandermonde)
    {
        if (!isVandermonde(A))
        {
            FatalErrorInFunction
                << "Source matrix not of Vandermonde type." << nl
                << abort(FatalError);
        }
    }

    for (label i = 0; i < n_; i++)
    {
        (*this)[i] = A[1][i];
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Vandermonde::~Vandermonde()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::Vandermonde::isVandermonde(const scalarSquareMatrix& A) const
{
    for (label i = 0; i < n_; i++)
    {
        for (label j = 0; j < n_; j++)
        {
            if (A[i][j] != pow(A[1][j], i))
            {
                return false;
            }
        }
    }

    return true;
}

void Foam::Vandermonde::solve
(
    scalarDiagonalMatrix& x,
    const scalarDiagonalMatrix& source
)
{
    if (n_ == 1)
    {
        x[0] = source[0];
        return;
    }

    scalarDiagonalMatrix c(n_, 0.0);

    // Calculate coefficients
    c[n_ - 1] = -(*this)[0];

    for (label i = 1; i < n_; i++)
    {
        const scalar xi = -(*this)[i];

        for (label j = n_ - i - 1; j < n_ - 1; j++)
        {
            c[j] += xi*c[j + 1];
        }

        c[n_ - 1] += xi;
    }

    // Solve system
    for (label i = 0; i < n_; i++)
    {
        const scalar xi = (*this)[i];

        scalar t = 1.0;
        scalar r = 1.0;
        scalar s = source[n_ - 1];

        for (label j = n_ - 1; j > 0; j--)
        {
            r = c[j] + r*xi;
            s += r*source[j - 1];
            t = r + t*xi;
        }

        x[i] = s / t;
    }
}

Foam::scalarSquareMatrix Foam::Vandermonde::invert()
{
    scalarSquareMatrix inverse(n_);
    scalarDiagonalMatrix source(n_, 0.0);
    scalarDiagonalMatrix x(n_);

    for (label i = 0; i < n_; i++)
    {
        // Build source vector
        if (i > 0)
        {
            source[i-1] = 0.0; 
        }
        
        source[i] = 1.0;

        // Solve
        solve(x, source);

        // Copy solution column
        for (label j = 0; j < n_; j++)
        {
            inverse[j][i] = x[j];
        }
    }

    return inverse;
}


// ************************************************************************* //
