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
    Copyright (C) 2019-2020 Alberto Passalacqua
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
        for (label i = 0; i < n_; i++)
        {
            for (label j = 0; i < n_; j++)
            {
                if (A[i][j] != pow(A[1][j], i))
                {
                    FatalErrorInFunction
                        << "Source matrix not of Vandermonde type." << nl
                        << abort(FatalError);
                }
            }
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

void Foam::Vandermonde::solve
(
    scalarDiagonalMatrix& x,
    const scalarDiagonalMatrix& source
)
{
    scalarDiagonalMatrix c(n_, 0.0);

    scalar t = 1.0;
    scalar r = 1.0;
    scalar s = 0.0;
    scalar xx = 0.0;

    if (n_ == 1)
    {
        x[0] = source[0];
    }
    else
    {
        c[n_ - 1] = -(*this)[0];

        for (label i = 1; i < n_; i++)
        {
            xx = -(*this)[i];

            for (label j = n_ - i - 1; j < n_ - 1; j++)
            {
                c[j] += xx*c[j + 1];
            }
            
            c[n_ - 1] += xx;
        }

        for (label i = 0; i < n_; i++)
        {
            xx = (*this)[i];
            t = 1.0;
            r = 1.0;
            s = source[n_ - 1];

            for (label j = n_ - 1; j > 0; j--)
            {
                r = c[j] + r*xx;
                s += r*source[j - 1];
                t = r + t*xx;
            }

            x[i] = s/t;
        }
    }
}

Foam::scalarSquareMatrix Foam::Vandermonde::inv()
{
    scalarSquareMatrix inv(n_);
    scalarDiagonalMatrix source(n_);
    scalarDiagonalMatrix x(n_);

    for (label i = 0; i < n_; i++)
    {
        for (label j = 0; j < n_; j++)
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

        for (label j = 0; j < n_; j++)
        {
            inv[j][i] = x[j];
        }
    }

    return inv;
}


// ************************************************************************* //
