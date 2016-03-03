/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 Alberto Passalacqua
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
    Foam::univariateMomentSet

Description
    Stores and operates on a set of univariate moments.
    
    References
    \verbatim
        "Modified moments and Gaussian quadratures"
        J C Wheeler
        Rocky Mountain Journal of Mathematics
        Volume 4, Pages 287-296, 1974
    \endverbatim
    
    \verbatim
        "Numerical Recipes in C++"
        W H Press, S A Teukolsky, W Vetterling, B P Flannery
        Cambridge University Press, 2007
    \endverbatim
    
    \verbatim
        "The theory of canonical moments with applications in Statistics, 
        Probability and Analysis"
        H Dette, W. J. Studden 
        Wiley & Sons, 1997
    \endverbatim

SourceFiles
    univariateMomentSet.C
    univariateMomentSetI.H

\*---------------------------------------------------------------------------*/

#ifndef univariateMomentSet_H
#define univariateMomentSet_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class univariateMomentSet Declaration
\*---------------------------------------------------------------------------*/

class univariateMomentSet
: 
    public scalarDiagonalMatrix
{
    // Private data

        //- Number of moments
        label nMoments_;

        //- Position of fist negative zeta
        label negativeZeta_;

        //- Bool to indicate if moments were inverted
        bool inverted_;

        //- Bool to indicate if the full moment set is realizable
        //  It is set to false also if a subset of moments is realizable
        bool realizable_;

        //- Bool to indicate moment realizability was checked
        bool realizabilityChecked_;
        
        //- Bool to indicate quadrature setup was completed
        bool quadratureSetUp_;

        //- Number of invertible moments
        label nInvertibleMoments_;
        
        //- Number of realizable moments
        label nRealizableMoments_;
        
        //- Number of quadrature nodes
        label nNodes_;
        
        //- Quadrature weights
        scalarDiagonalMatrix weights_;
        
        //- Quadrature abscissae
        scalarDiagonalMatrix abscissae_;
        
    // Private member functions
        
        //- Check moment realizability
        void checkRealizability();
        
        //- Setup quadrature storage based on realizability
        void setupQuadrature(bool clear = false);

public:

    // Constructors

        //- Construct from number of moments and initial value
        univariateMomentSet
        (
            const label nMoments,
            const scalar initValue
        );

        //- Construct from diagonalMatrix
        univariateMomentSet
        (
            const scalarDiagonalMatrix& m
        );

    //- Destructor
    virtual ~univariateMomentSet();


    // Member Functions

        //- Return quadrature abscissae
        inline const scalarDiagonalMatrix& abscissae() const;
        
        //- Invert the set of moments to compute weights and abscissae
        void invert();

        //- Update moments from weights and abscissae
        void update();

        //- Returns the status of inversion
        inline bool isInverted() const;
        
        //- Return realizability status of the moment set
        inline bool isRealizable();

        //- Return the number of invertible moments
        inline label nInvertibleMoments();

        //- Return the number of realizable moments
        inline label nRealizableMoments();

        //- Return the number of moments
        inline label nMoments() const;
        
        //- Return the number of nodes
        inline label nNodes() const;
        
        //- Return quadrature weights
        inline const scalarDiagonalMatrix& weights() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
