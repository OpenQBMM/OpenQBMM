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
    // Private Member Functions

        //- Check if a square matrix is Vandermonde
        bool isVandermonde(const scalarSquareMatrix& A) const;

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


    // Member Functions

        //- Solve the Vandermonde linear system with the given source vector x
        void solve
        (
            scalarDiagonalMatrix& x,
            const scalarDiagonalMatrix& source
        ) const;

        //- Invert Vandermonde matrix
        scalarSquareMatrix invert() const;


    // Static Member Functions

        //- The conditional solve of the quadratures that carry a quantity
        //  on the nodes of a univariate one: the matrix that takes the
        //  moments of orders 0 to m - 1 of the quantity, weighted by the
        //  nodes, to its value at each node,
        //
        //      diag(1/w) V^-1,
        //
        //  over the m nodes whose weight is above systemWeight, whose
        //  positions among the abscissae are returned in nodes.
        //
        //  A node below systemWeight is left out of the system, which is
        //  exact only for a node that carries nothing: one that carries a
        //  little still contributes to the moments, and leaving it out
        //  would give its share to the others. So systemWeight is zero
        //  where the moments are those of the whole quadrature. A node kept
        //  whose weight is not above solvedWeight is solved for, and so
        //  kept apart from the others, but its value is not trusted and its
        //  row is zero.
        //
        //  Nodes that carry nothing are left out because two of them at the
        //  same abscissa made the system singular. Returns false when two
        //  of the nodes kept share one, which a quadrature that has
        //  collapsed can leave.
        static bool conditionalInverse
        (
            const UList<scalar>& abscissae,
            const UList<scalar>& weights,
            const scalar systemWeight,
            const scalar solvedWeight,
            labelList& nodes,
            scalarSquareMatrix& invVR
        );


    // Member Operators

        //- Return the Vandermonde matrix
        inline scalarSquareMatrix operator()() const;

        //- Return a given element of the Vandermonde matrix
        inline scalar operator()(label i, label j) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VandermondeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
