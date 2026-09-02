/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
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
    Foam::univariateMomentSet

Description
    Stores a vector of integer moments of a univariate distribution, and
    provides the realizability check for the moment vector.

    The realizability check is the Wheeler algorithm, which builds the
    coefficients alpha_k and beta_k of the recurrence relationship of the
    monic orthogonal polynomials associated to the moment vector, walked
    together with the chain of the zeta_k values

    \f[
        \zeta_0 = \alpha_0, \quad
        \zeta_{2j-1} = \frac{\beta_j}{\zeta_{2j-2}}, \quad
        \zeta_{2j} = \alpha_j - \zeta_{2j-1}
    \f]

    which is the quantity whose positivity is equivalent to realizability for
    measures with support on R+. The walk stops at the first non-positive
    value, so the check also returns how many moments are realizable.

    The support of the measure selects the test applied at each step:

      - R: only the beta_k are required to be positive
      - R+: every zeta_k is required to be positive
      - [0, 1]: the canonical moments derived from the zeta_k are required to
        belong to [0, 1]

SourceFiles
    univariateMomentSet.C
    univariateMomentSetRealizability.C
    univariateMomentSetI.H

\*---------------------------------------------------------------------------*/

#ifndef univariateMomentSet_H
#define univariateMomentSet_H

#include "momentSet.H"
#include "scalarMatrices.H"
#include "supportType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class univariateMomentSet Declaration
\*---------------------------------------------------------------------------*/

class univariateMomentSet : public momentSet
{
public:

    // Public data types

        //- Outcome of the realizability check.
        //  Grouping the outcome makes every exit of the check assign the
        //  complete state, instead of leaving part of it stale.
        struct realizabilityStatus
        {
            //- Number of realizable moments
            label nRealizableMoments = 0;

            //- Bool to indicate if the full moment set is realizable
            bool fullyRealizable = true;

            //- Bool to indicate if a subset of the full moment set is
            //  realizable
            bool subsetRealizable = true;

            //- Bool to indicate if the moment set is on the boundary of the
            //  moment space
            bool onMomentSpaceBoundary = false;

            //- Bool to indicate the realizability check was performed
            bool checked = false;
        };


private:

    // Private data

        //- Number of additional quadrature points for the recurrence relation
        label nAdditionalQuadraturePoints_;

        //- Alpha coefficients of the recurrence relation
        scalarList alpha_;

        //- Beta coefficients of the recurrence relation
        scalarList beta_;

        //- Zeta values to check moment realizability
        scalarList zeta_;

        //- Canonical moments associated to the moment vector
        scalarList canonicalMoments_;

        //- Matrix used to build the recurrence relation
        scalarRectangularMatrix zRecurrence_;

        //- Outcome of the realizability check
        realizabilityStatus status_;


     // Private member functions

        //- Check moment realizability
        void checkRealizability
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Manage the moment vectors whose realizability is established
        //  without building the recurrence relationship: a null or negative
        //  zero-order moment, and a moment vector with only two moments.
        //  Returns true if the realizability of the moment set was settled.
        bool checkTrivialCases(bool fatalErrorOnFailedRealizabilityTest);

        //- Set the first row of the Wheeler table from the normalized
        //  moments, and the first coefficients of the recurrence
        //  relationship
        void initialiseRecurrence();

        //- Check zeta_0, which carries its own handling of the degenerate
        //  moment vector. Returns true if the realizability of the moment
        //  set was settled.
        bool checkFirstZeta(bool fatalErrorOnFailedRealizabilityTest);

        //- Walk the recurrence relationship and the chain of the zeta_k with
        //  the Wheeler algorithm, stopping at the first value that violates
        //  realizability. Returns the position of that value, or the number
        //  of zeta_k determined by the moments if the walk completed.
        label buildRecurrence();

        //- Set the realizability status from the position at which the walk
        //  of the recurrence relationship stopped
        void classify(const label zetai, const bool complete);

        //- Compute the canonical moments for measures with support over
        //  [0, 1] from the first nZeta values of the zeta chain
        void checkCanonicalMoments(const scalarList& zeta, const label nZeta);

        //- Number of zeta_k determined by the moment vector
        inline label nZeta() const;

        //- Reset the realizability status and the quantities derived from
        //  the moments
        void resetStatus();

        //- Returns the number of elements of the alpha vector of the
        //  recurrence relation. This is the number of quadrature nodes
        //  obtained assuming all the moments are invertible, plus the number
        //  of additional quadrature nodes.
        label nAlphaRecurrence(const label nAdditionalQuadratureNodes) const;

        //- Returns the number of elements of the beta vector of the recurrence
        //  relation.
        label nBetaRecurrence(const label nAdditionalQuadratureNodes) const;

        //- Returns the number of elements of the zeta and canonical moment
        //  vectors. Room is reserved for the additional quadrature nodes,
        //  which extend the zeta chain beyond what the moments determine.
        label nZetaRecurrence(const label nAdditionalQuadratureNodes) const;


public:

    // Constructors

        //- Construct from number of moments and initial value
        univariateMomentSet
        (
            const label nMoments,
            const supportType& support,
            const scalar smallM0,
            const scalar smallZeta,
            const scalar initialValue = 0,
            const label nAdditionalQuadraturePoints = 0
        );

        //- Construct from scalarList
        univariateMomentSet
        (
            const scalarList& m,
            const supportType& support,
            const scalar smallM0,
            const scalar smallZeta,
            const label nAdditionalQuadraturePoints = 0
        );


    // Member Functions

        //- Returns alpha coefficients of the recurrence relation
        inline const scalarList& alphaRecurrence() const;

        //- Returns alpha coefficients of the recurrence relation
        inline scalarList& alphaRecurrence();

        //- Returns beta coefficients of the recurrence relation
        inline const scalarList& betaRecurrence() const;

        //- Returns beta coefficients of the recurrence relation
        inline scalarList& betaRecurrence();

        //- Returns the canonical moments if support is [0, 1]
        inline scalarList& canonicalMoments();

        //- Returns the canonical moments if support is [0, 1]
        inline const scalarList& canonicalMoments() const;

        //- Returns true if only the zero-order moment is non-null
        inline bool isDegenerate() const;

        //- Return true if the entire moment set is realizable
        inline bool isFullyRealizable
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Return true if some moments are realizable
        inline bool isRealizable
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Return true if a subset of moments is realizable
        inline bool isSubsetRealizable
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Return true if the moment set is on the boundary of the moment space
        inline bool isOnMomentSpaceBoundary
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Build the list of moment orders of a univariate moment set
        static labelListList makeUnivariateMomentOrders(const label nMoments);

        //- Const access to the moment of order momentOrder, without
        //  invalidating the realizability check
        inline const scalar& moment(const label momentOrder) const;

        //- Return the number of realizable moments
        inline label nRealizableMoments
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Update moments from the corresponding quadrature
        virtual void update
        (
            const scalarList& weights,
            const scalarList& abscissae
        );

        //- Update integer moments from the corresponding quadrature
        virtual void updateIntegerMoments
        (
            const scalarList& weights,
            const scalarList& abscissae
        );

        //- Set the alpha and beta coefficients of the recurrence relation
        //  from the first nZeta values of the zeta chain. Valid for measures
        //  with support over R+ and [0, 1], for which the zeta chain is the
        //  primary representation of the recurrence relationship.
        //  The coefficients are written into the lists of the caller, and
        //  those of the moment set are left untouched, so that a quadrature
        //  can extend the recurrence relationship without altering the one
        //  the moment vector determines.
        void zetasToRecurrence
        (
            const label nZeta,
            scalarList& alpha,
            scalarList& beta
        ) const;

        //- Const access to element
        inline const scalar& operator[](const label index) const;

        //- Non-const access to element, reinitializing realizaibility check
        inline scalar& operator[](const label index);

        //- Const access to the moment of order momentOrder
        inline const scalar& operator()(const label momentOrder) const;

        //- Non-const access to the moment of order momentOrder
        inline scalar& operator()(const label momentOrder);

        //- Resizes moment set
        void resize(const label newSize);

        //- Set size of moment set
        void setSize(const label newSize);

        //- Return support
        inline const supportType& support() const;

        //- Returns zeta_k values for the moment set
        inline scalarList& zetas();

        //- Returns zeta_k values for the moment set
        inline const scalarList& zetas() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
