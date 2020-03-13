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

Class
    Foam::Vandermonde

Description
    Stores the second row of a Vandermonde matrix, and solves the associated
    linear system.

    References
    \verbatim
        William H. Press, Saul A. Teukolsky,
        William T. Vetterling, and Brian P. Flannery. 1992.
        "Numerical Recipes in C (2nd Ed.): The Art of Scientific Computing."
        Cambridge University Press, New York, NY, USA.
    \endverbatim

SourceFiles
    Vandermonde.C

\*---------------------------------------------------------------------------*/

#ifndef Vandermonde_H
#define Vandermonde_H

#include "scalarMatrices.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class Vandermonde Declaration
\*---------------------------------------------------------------------------*/

class Vandermonde
:
    public scalarDiagonalMatrix
{
    // Private data

        //- Number of rows and columns (square matrix)
        const label n_;

public:

    // Constructors

        //- Construct from a given vector containing the second row of the 
        //  Vandermonde matrix
        Vandermonde
        (
            const scalarDiagonalMatrix& A
        );

        //- Construct from a square Vandermonde matrix
        Vandermonde
        (
            const scalarSquareMatrix& A,
            const bool checkVandermonde = false
        );


    //- Destructor
    ~Vandermonde();


    // Member Functions

        //- Solve the Vandermonde linear system with the given source vector x
        void solve
        (
            scalarDiagonalMatrix& x,
            const scalarDiagonalMatrix& source
        );

        //- Invert Vandermonde matrix
        scalarSquareMatrix inv();


    // Member Operators

        //- Return the Vandermonde matrix
        inline scalarSquareMatrix operator()();

        //- Return a given element of the Vandermonde matrix
        inline scalar operator()(label i, label j);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VandermondeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
