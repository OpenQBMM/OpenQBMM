/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2026 Alberto Passalacqua
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
    Foam::generalizedMomentInversion

Description
    Generalized quadrature method of moments (GQMOM): a Gauss quadrature of
    more nodes than the moments determine, which reproduces the moments it
    is given and takes the rest of the recurrence relationship of its
    orthogonal polynomials from a model distribution.

    The moments determine the first coefficients of the recurrence
    relationship, a_i and b_i, or, for a support on R+ and on [0, 1], the
    first zeta_k, with a_0 = zeta_1, a_i = zeta_(2i) + zeta_(2i+1) and
    b_i = zeta_(2i-1) zeta_(2i), and the first canonical moments p_k, with
    zeta_k = p_k (1 - p_(k-1)). The generalized quadrature continues them
    with the law of the polynomials of a model distribution, scaled so that
    it starts from the last value the moments determine: a distribution of
    the model is reproduced exactly, its higher moments included, and a
    moment set on the boundary of the moment space, whose last value is
    zero, keeps the quadrature of the regular nodes.

    - R, Hermite polynomials of a Gaussian distribution, b_i = i sigma^2:
      a_i is the mean of the known ones, and
      \f[
          b_i = \left( \frac{i}{l} \right)^\nu b_l
      \f]
      with b_l the last known one and nu the exponent that controls the
      tails (nu = 1 gives the Gaussian ones, the default).

    - R+, ndfTypeRPlus gamma (the default), generalized Laguerre
      polynomials of x^alpha exp(-beta x), with alpha from the mean and the
      variance:
      \f[
          \zeta_{2k-1} = \frac{k + \alpha}{\beta}, \quad
          \zeta_{2k} = \frac{k}{\beta}
      \f]

    - R+, ndfTypeRPlus lognormal, Stieltjes-Wigert polynomials, with
      eta^2 = m_0 m_2/m_1^2:
      \f[
          \zeta_{2k-1} = e^\mu \eta^{4k-3}, \quad
          \zeta_{2k} = e^\mu \eta^{2k-1} (\eta^{2k} - 1)
      \f]

    - [0, 1], Jacobi polynomials of x^beta (1 - x)^alpha, with alpha and
      beta from p_1 and p_2:
      \f[
          p_{2k-1} = \frac{\beta + k}{2k + \alpha + \beta}, \quad
          p_{2k} = \frac{k}{2k + 1 + \alpha + \beta}
      \f]
      A canonical moment above the one of the model is scaled from one
      rather than from zero, so that it stays inside (0, 1).

    The values of odd index and those of even index are two chains, and each
    continues from its own last known value, with the ratio of the law
    between the index sought and the index of that value. How far each
    chain is known depends on the number of moments. With 2n + 1 moments,
    up to order 2n, which is the case the reference is written for, both
    end at the index n: zeta_(2n-1) and zeta_(2n), and b_n. With 2n
    moments, up to order 2n - 1, which is what a case that carries the
    moments of a Gauss quadrature of n nodes has, the even chain ends one
    index earlier, at zeta_(2n-2), as b does at b_(n-1), and is continued
    from there, zeta_(2n) included. Where only some of the moments are
    realizable, the chains end where the realizable ones do.

    Reference
    \verbatim
        "The generalized quadrature method of moments"
        R O Fox, F Laurent, A Passalacqua
        Journal of Aerosol Science
        Volume 167, 106096, 2023
    \endverbatim

SourceFiles
    generalizedMomentInversion.C
    generalizedMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef generalizedMomentInversion_H
#define generalizedMomentInversion_H

#include "univariateMomentInversion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class generalizedMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class generalizedMomentInversion
:
    public univariateMomentInversion
{
    // Private data

        // Exponent
        scalar nu_;

        // NDF type for R+ support
        word ndfTypeRPlus_;

        // Maximum number of nodes for generalized quadrature
        label nMaxNodes_;

        // Number of nodes from regular quadrature (internally calculated)
        label nRegularQuadratureNodes_;

        // Number of additional quadrature nodes for generalized quadrature
        label nAdditionalQuadratureNodes_;

        //- Number of realizable moments of the set being inverted, which
        //  says how many coefficients of the recurrence relationship the
        //  moments determine, and so where the closure starts
        label nRealizableMoments_;


    // Private member functions

        //- Returns true if the moment vector supports the extension of the
        //  recurrence relationship to the additional quadrature nodes
        bool canExtendRecurrence(const univariateMomentSet& moments) const;

        // Correct recurrence relation for R support
        void correctRecurrenceR
        (
            scalarList& alpha,
            scalarList& beta
        );

        // Correct recurrence relation for R+ support
        void correctRecurrenceRPlus
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta
        );

        // Correct recurrence relation for [0, 1] support
        void correctRecurrence01
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta
        );

public:

    //- Runtime type information
    TypeName("GQMOM");


    // Constructors

        //- Construct from univariateMomentSet
        generalizedMomentInversion
        (
            const dictionary& dict,
            const label nMaxNodes = 0
        );



    // Public member functions

        //- Modify the recurrence relation for fixed-point quadrature
        virtual void correctRecurrence
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta
        );

        //- Calculates the number of quadrature nodes
        virtual void calcNQuadratureNodes
        (
            univariateMomentSet& moments
        );

        //- Return the number of additional quadrature nodes. Generalized
        //  QMOM extends the recurrence relationship of the moments to the
        //  nodes the user asks for beyond the ones the moments determine
        virtual label nAdditionalQuadraturePoints(const label nMoments) const;

        //- Invert the set of moments to compute weights and abscissae
        void invert
        (
            univariateMomentSet& moments
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
