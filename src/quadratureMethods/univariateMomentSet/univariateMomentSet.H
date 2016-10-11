/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 Alberto Passalacqua
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
    Stores, checks realizability and calculates the quadrature representation of
    a measure from a set of univariate moments.

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

    \verbatim
        "Matrix measures, moment spaces and Favard’s theorem for the interval
        [0, 1] and [0,∞)"
        H Dette, W. J. Studden
        Linear algebra and its applications
        Issue 345, Pages 169-193, 2002
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

        //- Alpha coefficients of the recurrence relation
        scalarDiagonalMatrix alpha_;

        //- Beta coefficients of the recurrence relation
        scalarDiagonalMatrix beta_;

        //- Type of quadrature formula computed when inverting the moment set
        word quadratureType_;

        //- Type of support of the measure associated to the moment set
        word support_;

        //- Position of first negative zeta
        label negativeZeta_;

        //- Bool to indicate if moment set is degenerate, with only the
        //  zero-order moment realizable
        bool degenerate_;

        //- Bool to indicate if moments were inverted
        bool inverted_;

        //- Bool to indicate if the full moment set is realizable
        bool fullyRealizable_;

        //- Bool to indicate if a subset of the full moment set is realizable
        bool subsetRealizable_;

        //- Bool to indicate if the moment set is on the boundary of the
        //  moment space
        bool onMomentSpaceBoundary_;

        //- Bool to indicate moment realizability was checked
        bool realizabilityChecked_;

        //- Bool to indicate quadrature setup was completed
        bool quadratureSetUp_;

        //- Bool to force Gauss quadrature when using other quadrature
        //  formulae. Usually needed to preserve the largest possible
        //  number of moments.
        bool forceGauss_;

        //- Number of invertible moments
        label nInvertibleMoments_;

        //- Number of realizable moments
        label nRealizableMoments_;

        //- Number of quadrature nodes
        label nNodes_;

        //- Known abscissa for Gauss-Radau quadrature
        //  Optionally specified at construction. If not specified, it is
        //  assumed to be zero.
        scalar knownAbscissa_;

        //- Quadrature weights
        scalarDiagonalMatrix weights_;

        //- Quadrature abscissae
        scalarDiagonalMatrix abscissae_;

    // Private member functions

        //- Check moment realizability
        void checkRealizability();

        //- Checks canonical moments for measures with support over [0,1]
        void checkCanonicalMoments
        (
            const scalarDiagonalMatrix& zeta,
            const label nZeta
        );

        //- Computes the number of invertible moments
        inline void calcNInvertibleMoments();

        //- Setup quadrature storage based on realizability
        void setupQuadrature(bool clear = false, bool nullMomentSet = false);

public:

    // Constructors

        //- Construct from number of moments and initial value
        univariateMomentSet
        (
            const label nMoments,
            const scalar initValue,
            const word& quadratureType,
            const word& support,
            const scalar knownAbscissa = 0.0
        );

        //- Construct from diagonalMatrix
        univariateMomentSet
        (
            const scalarDiagonalMatrix& m,
            const word& quadratureType,
            const word& support,
            const scalar knownAbscissa = 0.0
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

        //- Returns degenerate status of the moment set
        inline bool isDegenerate();

        //- Returns the status of inversion
        inline bool isInverted() const;

        //- Return realizability status of the moment set
        inline bool isFullyRealizable();

        //- Return realizability status of the moment subsets
        inline bool isSubsetRealizable();

        //- Return true if the moment set is on the boundary of the moment space
        inline bool isOnMomentSpaceBoundary();

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

        //- Return support
        inline const word quadratureType() const;

        //- Return support
        inline const word support() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
