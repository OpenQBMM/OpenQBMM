/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
    Copyright (C) 2026 Alberto Passalacqua
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
    Foam::GolubWelsch

Description
    Gaussian quadrature from the coefficients of the three-term recurrence
    relationship of the monic orthogonal polynomials associated to a measure.

    The abscissae of the quadrature are the eigenvalues of the Jacobi matrix
    of the recurrence relationship

    \f[
        J =
        \begin{bmatrix}
            \alpha_0    & \sqrt{\beta_1} &                &        \\
            \sqrt{\beta_1} & \alpha_1    & \sqrt{\beta_2} &        \\
                        & \ddots         & \ddots         & \ddots \\
        \end{bmatrix}
    \f]

    and the weights are \f$ w_i = m_0 q_{0i}^2 \f$, where \f$ q_{0i} \f$ is
    the first component of the eigenvector associated to the i-th eigenvalue.

    The Jacobi matrix is symmetric and tridiagonal by construction, so it is
    stored as its diagonal and sub-diagonal, and the eigenvalue problem is
    solved by the QL algorithm with implicit shifts alone, without the
    Householder reduction a general symmetric matrix requires. Only the first
    row of the product of the Givens rotations is accumulated, because it is
    the only one the weights need.

    References
    \verbatim
        Golub, G. H., and Welsch, J. H. (1969).
        "Calculation of Gauss quadrature rules."
        Mathematics of Computation, 23(106), 221-230.
    \endverbatim

    The QL sweep follows the one of Foam::EigenMatrix, itself derived from the
    JAMA library, so that the two give the same answer on the same matrix.

SourceFiles
    GolubWelschI.H

\*---------------------------------------------------------------------------*/

#ifndef GolubWelsch_H
#define GolubWelsch_H

#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class GolubWelsch Declaration
\*---------------------------------------------------------------------------*/

class GolubWelsch
{
    // Private data

        //- Diagonal of the Jacobi matrix, overwritten with its eigenvalues
        scalarList d_;

        //- Sub-diagonal of the Jacobi matrix, consumed by the QL sweep
        scalarList e_;

        //- First row of the accumulated product of the Givens rotations
        scalarList q_;


    // Private member functions

        //- Grow the work lists to hold nNodes entries. They are only grown,
        //  never shrunk, so that repeated inversions of the same size, which
        //  is what a solver does, do not allocate.
        inline void resizeWork(const label nNodes);

        //- Eigenvalues of the symmetric tridiagonal matrix stored in d_ and
        //  e_, by the QL algorithm with implicit shifts, accumulating the
        //  rotations on q_ alone
        inline void symmTridiagonalQL(const label nNodes);


public:

    // Constructors

        //- Construct null
        inline GolubWelsch();

        //- Construct with room for nNodes
        inline explicit GolubWelsch(const label nNodes);


    // Member Functions

        //- Compute the Gaussian quadrature of nNodes nodes determined by the
        //  coefficients of the recurrence relationship and by the moment of
        //  order zero. Only beta[1] to beta[nNodes - 1] are read: beta[0] is
        //  the moment of order zero of the measure and does not enter the
        //  Jacobi matrix.
        inline void quadrature
        (
            const UList<scalar>& alpha,
            const UList<scalar>& beta,
            const label nNodes,
            const scalar m0,
            UList<scalar>& abscissae,
            UList<scalar>& weights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GolubWelschI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
