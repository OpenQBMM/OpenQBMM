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

        //- Dimesions of the matrix
        const label m_;

public:

    // Constructors

        //- Construct from a given vector
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

        //- Solve the Vandermonde linear system with the given source vector
        void solve
        (
            scalarDiagonalMatrix& x,
            const scalarDiagonalMatrix& source
        );

        //- Invert Vandermonde matrix
        scalarSquareMatrix inv();


    // Member Operators

        //- Return the complete Vandermonde matrix
        inline scalarSquareMatrix operator()();

        //- Return a given component of the Vandermonde matrix
        inline scalar operator()(label i, label j);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VandermondeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
