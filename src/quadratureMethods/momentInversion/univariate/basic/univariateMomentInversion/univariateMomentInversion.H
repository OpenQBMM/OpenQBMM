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
    Foam::univariateMomentInversion

Description
    Abstract class for univariate quadrature on a moment set. The quadrature
    approach is selected at runtime.

    Every inversion reads from its dictionary

    \table
        Property          | Description                          | Default
        smallM0           | Smallest zero-order moment inverted  | SMALL
        smallZeta         | Smallest accepted zeta_k             | 0
        boundaryTolerance | On [0, 1], how close to zero or to \
                            one a canonical moment is taken to \
                            be at it                             | 1e-12
    \endtable

    A moment set on [0, 1] with a canonical moment at an end is on the
    boundary of the moment space, where its measure is unique, and is given
    that measure whatever quadrature was asked for (see invertOnBoundary).
    The default of boundaryTolerance is what round-off leaves of a canonical
    moment that is exactly zero or one, as the one of two pure streams is.
    A larger value also treats as pure a stream that is almost so: a node
    the moments put within about that distance of an end of the interval is
    moved to the end. The moments of order zero and one are kept in either
    case.

    A canonical moment beyond that band, below zero or above one, belongs to
    no measure: the moments of that order and above are then replaced by
    those of the measure of the ones below it, which is the most the
    inversion can return. Whether the last set inverted was such a set is
    kept in boundaryBeyond, for a caller that inverts a field to count.

SourceFiles
    univariateMomentInversion.C
    univariateMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef univariateMomentInversion_H
#define univariateMomentInversion_H

#include "scalarList.H"
#include "scalarMatrices.H"
#include "runTimeSelectionTables.H"
#include "univariateMomentSet.H"
#include "GolubWelsch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class univariateMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class univariateMomentInversion
{
    // Private data

        //- Smallest accepted value for m0
        scalar smallM0_;

        //- Smallest accepted value for zeta_k
        scalar smallZeta_;

        //- On [0, 1], how close to an end of the interval a canonical
        //  moment is taken to be at it, so that the moment set is on the
        //  boundary of the moment space
        scalar boundaryTolerance_;

        //- Whether the moment set of the last inversion had a canonical
        //  moment beyond an end of [0, 1]
        bool boundaryBeyond_;

        //- How far beyond that end it lay
        scalar boundaryExcess_;


protected:

    // Protected data

        //- Number of invertible moments for the user-selected quadrature
        label nInvertibleMoments_;

        //- Number of quadrature nodes
        label nNodes_;

        //- Quadrature abscissae
        scalarList abscissae_;

        //- Quadrature weights
        scalarList weights_;

        //- Gaussian quadrature of the recurrence relationship. It holds the
        //  work arrays of the eigenvalue problem, so that repeated
        //  inversions of the same size do not allocate.
        GolubWelsch golubWelsch_;

        //- Work copy of the alpha coefficients of the recurrence relation.
        //  correctRecurrence builds the coefficients of the quadrature, which
        //  are not the ones the moment vector determines: Gauss-Radau and
        //  Gauss-Lobatto correct the last of them to fix a known abscissa,
        //  and GQMOM extends them to the additional nodes. The coefficients
        //  of the moment set are left untouched, so the work is done on a
        //  copy, kept as a member to avoid allocating it on every inversion.
        //  Note that GQMOM does extend the zeta chain of the moment set in
        //  place, which the realizability check rebuilds for the next cell.
        scalarList alpha_;

        //- Work copy of the beta coefficients of the recurrence relation
        scalarList beta_;


    // Protected member functions

        //- Build the coefficients of the recurrence relationship of the
        //  quadrature, in the alpha_ and beta_ work lists
        void correctedRecurrence
        (
            univariateMomentSet& moments
        );


        //- The measure of a moment set on the boundary of the moment space
        //  of [0, 1], which is the only one with those moments, so that no
        //  rule has anything to choose: returns false for a set in the
        //  interior, or on another support.
        //
        //  The zeta_k of such a measure stop: with the first canonical
        //  moment at an end being p_N, zeta_N is zero where p_N is zero,
        //  and zeta_(N+1) is where p_N is one, zeta_N being then 1 -
        //  p_(N-1). The Jacobi matrix the zeta_k give, alpha_j = zeta_2j +
        //  zeta_(2j+1) and beta_j = zeta_(2j-1) zeta_2j, stops with them,
        //  and its eigenvalues are the nodes: those of a Gauss rule, or a
        //  node at zero, at one or at both with the interior ones,
        //  whichever the moments say. Two pure streams, m_k = m_1 for every
        //  k, have p_2 = 1 and come back as nodes at zero and one of
        //  weights 1 - m_1 and m_1; a Gauss rule gives them one node at
        //  the mean, which is the two streams mixed.
        bool invertOnBoundary(univariateMomentSet& moments);

        //- Calculate weights and abscissae
        void calcQuadrature(const univariateMomentSet& moments);

        //- Evaluate the monic orthogonal polynomials of the recurrence
        //  relationship at x, up to degree nNodes_ - 1, and return the value
        //  of the polynomial of degree nNodes_ - 2 in pMinus1. Used by
        //  Gauss-Radau and Gauss-Lobatto to fix a known abscissa.
        scalar orthogonalPolynomial
        (
            const scalarList& alpha,
            const scalarList& beta,
            const scalar x,
            scalar& pMinus1
        ) const;


public:

    //- Runtime type information
    TypeName("univariateMomentInversion");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMomentInversion,
            dictionary,
            (
                const dictionary& dict,
                const label nMaxNodes
            ),
            (dict, nMaxNodes)
        );


    // Constructors

        //- Construct from univariateMomentSet
        univariateMomentInversion
        (
            const dictionary& dict,
            const label nMaxNodes = 0
        );

        //- Disallow default bitwise copy construct
        univariateMomentInversion(const univariateMomentInversion&) = delete;


    //- Destructor
    virtual ~univariateMomentInversion() = default;


    // Selectors

        //- The number of quadrature nodes the inversion the dictionary
        //  selects builds from nMoments moments, which is the number of
        //  nodes a case of that many moments carries
        static label nQuadratureNodes
        (
            const dictionary& dict,
            const label nMoments
        );

        static autoPtr<univariateMomentInversion> New
        (
            const dictionary& dict,
            const label nMaxNodes = 0
        );


    // Member Functions

        //- Calculates the number of quadrature nodes
        virtual void calcNQuadratureNodes
        (
            univariateMomentSet& moments
        ) = 0;

        //- Modify the recurrence relation for fixed-point quadrature
        virtual void correctRecurrence
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta
        ) = 0;

        //- Invert the set of moments to compute weights and abscissae
        virtual void invert
        (
            univariateMomentSet& moments
        );

        //- Return the number of quadrature nodes, in excess of the ones a
        //  Gauss quadrature builds from nMoments moments, for which the
        //  moment set has to reserve room in the coefficients of the
        //  recurrence relationship and in the chain of the zeta_k
        virtual label nAdditionalQuadraturePoints(const label nMoments) const;

        //- Return quadrature abscissae
        inline const scalarList& abscissae() const;

        //- Return the number of nodes
        inline label nNodes() const;

        //- Return quadrature weigths
        inline const scalarList& weights() const;

        //- Return the smallest accepted value for m0
        inline scalar smallM0() const;

        //- Return the smallest accepted value for zeta_k
        inline scalar smallZeta() const;

        //- Return how close to an end of [0, 1] a canonical moment is
        //  taken to be at it
        scalar boundaryTolerance() const
        {
            return boundaryTolerance_;
        }

        //- Whether the last moment set inverted had a canonical moment
        //  beyond an end of [0, 1], and so was the moment set of no
        //  measure. What was returned is the measure of the moments below
        //  that one. A caller that inverts a field counts these and
        //  reports them.
        bool boundaryBeyond() const
        {
            return boundaryBeyond_;
        }

        //- How far beyond an end of [0, 1] the canonical moment of the last
        //  moment set lay: what the report of a field says, so that the
        //  round-off the moments carry out of their own equations is told
        //  apart from a moment set that is wrong
        scalar boundaryExcess() const
        {
            return boundaryExcess_;
        }

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const univariateMomentInversion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentInversionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
