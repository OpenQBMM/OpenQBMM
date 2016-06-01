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


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Vandermonde::~Vandermonde()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::Vandermonde::solve
(
    const scalarDiagonalMatrix& b,
    scalarDiagonalMatrix& x
)
{
    const scalarDiagonalMatrix& A = *this;
    label n = this->size();
    scalarDiagonalMatrix c(n, scalar(0));
    scalar t = scalar(1);
    scalar r = scalar(1);
    scalar s = scalar(0);
    scalar xx = scalar(0);
    
    if (n == 1)
    {
        x[0] = b[0];
    }
    
    else
    {
        c[n-1] = -A[0];
        
        for (label i = 1; i < n; i++)
        {
            xx = -A[i];
            
            for (label j = (n-i-1); j < (n-1); j++)
            {
                c[j] += xx*c[j+1];
            }
            c[n-1] += xx;
        }
        
        for (label i = 0; i < n; i++)
        {
            xx = A[i];
            t = scalar(1);
            r = scalar(1);
            s = b[n-1];
            
            for (label j = (n-1); j > 0; j--)
            {
                r = c[j] + xx*r;
                s += b[j-1]*r;
                t = xx*t + r;
            }
            x[i] = s/t;
        }
    }
    return;
}

void Foam::Vandermonde::invert(scalarSquareMatrix& invA)
{    
    label n = this->size();
    scalarSquareMatrix identity(n, scalar(0));
    
    scalarDiagonalMatrix b(n);
    scalarDiagonalMatrix x(n);
    
    for (label i = 0; i < n; i++)
    {
        identity(i,i) = scalar(1);
        
        for (label j = 0; j < n; j++)
        {
            b[j] = identity(j,i);
        }
        solve(b, x);

        for (label j = 0; j < n; j++)
        {
            invA(j,i) = x[j];
        }
    }
}

// ************************************************************************* //
