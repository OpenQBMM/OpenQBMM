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
    const scalarSquareMatrix& A,
    const bool checkVandermonde
)
:
    scalarDiagonalMatrix(A.m())
{
    if (checkVandermonde)
    {
        if (!isVandermonde(A))
        {
            FatalErrorInFunction
                << "Source matrix not of Vandermonde type." << nl
                << exit(FatalError);
        }
    }

    const label n = this->size();

    if (n > 1)
    {
        for (label i = 0; i < n; i++)
        {
            (*this)[i] = A[1][i];
        }
    }
    else if (n == 1)
    {
        // A 1x1 Vandermonde matrix is [1] regardless of the base value,
        // so there is no second row to read it from; solve()/invert()/
        // operator() never consult it in this case either.
        (*this)[0] = 0.0;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::Vandermonde::isVandermonde(const scalarSquareMatrix& A) const
{
    const label n = this->size();

    for (label j = 0; j < n; j++)
    {
        const scalar base = (n > 1) ? A[1][j] : 0.0;
        scalar expectedPower = 1.0;

        for (label i = 0; i < n; i++)
        {
            // Relative tolerance: expectedPower grows like base^i, so a
            // fixed absolute tolerance rejects perfectly valid matrices
            // as soon as the abscissae move away from the unit scale.
            const scalar tol = ROOTSMALL*max(mag(expectedPower), 1.0);

            if (mag(A[i][j] - expectedPower) > tol)
            {
                return false;
            }

            expectedPower *= base;
        }
    }

    return true;
}

void Foam::Vandermonde::solve
(
    scalarDiagonalMatrix& x,
    const scalarDiagonalMatrix& source
) const
{
    const label n = this->size();

    if (source.size() != n)
    {
        FatalErrorInFunction
            << "Source vector size (" << source.size()
            << ") does not match matrix size (" << n << ")" << nl
            << exit(FatalError);
    }

    if (x.size() != n)
    {
        FatalErrorInFunction
            << "Solution vector size (" << x.size()
            << ") does not match matrix size (" << n << ")" << nl
            << exit(FatalError);
    }

    if (n == 1)
    {
        x[0] = source[0];
        return;
    }

    // A repeated (or nearly repeated) abscissa makes the Vandermonde
    // matrix singular. Detect it directly here: the natural scale of the
    // accumulated product t below is (abscissa)^(n - 1), which is not
    // comparable to a fixed absolute tolerance such as VSMALL.
    for (label i = 0; i < n; i++)
    {
        for (label j = i + 1; j < n; j++)
        {
            const scalar diff = mag((*this)[i] - (*this)[j]);
            const scalar scale = max(mag((*this)[i]), mag((*this)[j]));

            if (diff < ROOTSMALL*max(scale, 1.0))
            {
                FatalErrorInFunction
                    << "Near-singular Vandermonde matrix: abscissae " << i
                    << " and " << j << " are not sufficiently distinct."
                    << nl
                    << "  Value " << i << ": " << (*this)[i] << nl
                    << "  Value " << j << ": " << (*this)[j] << nl
                    << exit(FatalError);
            }
        }
    }

    scalarDiagonalMatrix c(n, 0.0);

    // Calculate coefficients
    c[n - 1] = -(*this)[0];

    for (label i = 1; i < n; i++)
    {
        const scalar xi = -(*this)[i];

        for (label j = n - i - 1; j < n - 1; j++)
        {
            c[j] += xi*c[j + 1];
        }

        c[n - 1] += xi;
    }

    // Solve system
    for (label i = 0; i < n; i++)
    {
        const scalar xi = (*this)[i];

        scalar t = 1.0;
        scalar r = 1.0;
        scalar s = source[n - 1];

        for (label j = n - 1; j > 0; j--)
        {
            r = c[j] + r*xi;
            s += r*source[j - 1];
            t = r + t*xi;
        }

        x[i] = s / t;
    }
}

Foam::scalarSquareMatrix Foam::Vandermonde::invert() const
{
    const label n = this->size();

    scalarSquareMatrix inverse(n);
    scalarDiagonalMatrix source(n, 0.0);
    scalarDiagonalMatrix x(n);

    for (label i = 0; i < n; i++)
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
        for (label j = 0; j < n; j++)
        {
            inverse[j][i] = x[j];
        }
    }

    return inverse;
}


// ************************************************************************* //
