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
    Stores the second row of a Vandermonde form. Solves a system of equations
    based on a square Vandermonde form with powers from 0 to n-1.

    References
    \verbatim
        William H. Press, Saul A. Teukolsky, 
        William T. Vetterling, and Brian P. Flannery. 1992. 
        "Numerical Recipes in C (2nd Ed.): The Art of Scientific Computing."
        Cambridge University Press, New York, NY, USA.

    \endverbatim
    
SourceFiles
    VandermondeI.H
    Vandermonde.C

\*---------------------------------------------------------------------------*/

#ifndef Vandermonde_H
#define Vandermonde_H

#include "PtrList.H"
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

// Public Data
public:

    // Constructors

        //- Construct for a given vector
        Vandermonde
        (
            const scalarDiagonalMatrix& A
        );


    //- Destructor
    ~Vandermonde();


    // Member Functions

        //- Solve Ax=b where A is the List of Vandermonde components
        void solve
        (
            const scalarDiagonalMatrix& b,
            scalarDiagonalMatrix& x
        );
            
        //- Invert Vancermonde system
        void invert(scalarSquareMatrix& invA);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
