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

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::scalarSquareMatrix Foam::Vandermonde::operator()()
{
    scalarSquareMatrix tmpV(n_);

    for (label j = 0; j < n_; j++)
    {
        // Cache base value
        const scalar base = (*this)[j];

        // Initialize iterative calculation of integer power
        scalar power = 1.0;
        tmpV[0][j] = power;
        
        // Iteratively compute matrix elements
        for (label i = 1; i < n_; i++)
        {
            power *= base;
            tmpV[i][j] = power;
        }
    }

    return tmpV;
}

Foam::scalar Foam::Vandermonde::operator()(label i, label j)
{
    #ifdef FULLDEBUG
    if (i >= n_ || j >= n_ || i < 0 || j < 0)
    {
        FatalErrorInFunction
            << "Index out of bounds: (" << i << "," << j 
            << ") for matrix size " << n_ << "x" << n_ << nl
            << abort(FatalError);
    }
    #endif
    
    // Handle trivial cases
    if (i == 0) 
    {
        return 1.0;
    }

    if (i == 1) 
    {
        return (*this)[j];
    }

    // Compute element using iterative calculation of integer power
    scalar result = (*this)[j];

    for (label k = 1; k < i; k++)
    {
        result *= (*this)[j];
    }
    
    return result;
}


// ************************************************************************* //
