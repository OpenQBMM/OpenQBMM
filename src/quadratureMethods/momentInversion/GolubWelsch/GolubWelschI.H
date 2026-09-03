/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::GolubWelsch::GolubWelsch()
:
    d_(),
    e_(),
    q_()
{}


inline Foam::GolubWelsch::GolubWelsch(const label nNodes)
:
    d_(nNodes, Zero),
    e_(nNodes, Zero),
    q_(nNodes, Zero)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::GolubWelsch::resizeWork(const label nNodes)
{
    if (d_.size() < nNodes)
    {
        d_.setSize(nNodes);
        e_.setSize(nNodes);
        q_.setSize(nNodes);
    }
}


inline void Foam::GolubWelsch::symmTridiagonalQL(const label nNodes)
{
    // A QL sweep converges in a handful of iterations per eigenvalue, so this
    // bound is only reached if the matrix is not the one the algorithm
    // assumes. It is here because the sweep has no other way out.
    const label maxIterations = 50;

    const scalar eps = pow(2.0, -52.0);

    scalar f = 0.0;
    scalar tst1 = 0.0;

    for (label l = 0; l < nNodes; l++)
    {
        // Look for a small sub-diagonal element
        tst1 = Foam::max(tst1, mag(d_[l]) + mag(e_[l]));

        label m = l;

        while (m < nNodes)
        {
            if (mag(e_[m]) <= eps*tst1)
            {
                break;
            }

            m++;
        }

        // If m == l, d_[l] is an eigenvalue, otherwise iterate. The sweep
        // reads d_[l + 1] and e_[l + 1], which m > l keeps in range
        if (m > l)
        {
            label iteration = 0;

            do
            {
                if (++iteration > maxIterations)
                {
                    FatalErrorInFunction
                        << "The QL algorithm did not converge to the "
                        << "eigenvalues of the Jacobi matrix." << nl
                        << "    Number of quadrature nodes: " << nNodes << nl
                        << "    Eigenvalue: " << l << nl
                        << "    Iterations: " << maxIterations << nl
                        << "    Diagonal: " << d_ << nl
                        << "    Sub-diagonal: " << e_ << nl
                        << exit(FatalError);
                }

                // Compute the implicit shift
                scalar g = d_[l];
                scalar p = (d_[l + 1] - g)/(2.0*e_[l]);
                scalar r = Foam::hypot(p, scalar(1));

                if (p < 0)
                {
                    r = -r;
                }

                d_[l] = e_[l]/(p + r);
                d_[l + 1] = e_[l]*(p + r);

                const scalar dl1 = d_[l + 1];

                scalar h = g - d_[l];

                for (label i = l + 2; i < nNodes; i++)
                {
                    d_[i] -= h;
                }

                f += h;

                // Implicit QL transformation
                p = d_[m];

                scalar c = 1.0;
                scalar c2 = c;
                scalar c3 = c;
                scalar s = 0.0;
                scalar s2 = 0.0;

                const scalar el1 = e_[l + 1];

                for (label i = m - 1; i >= l; i--)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c*e_[i];
                    h = c*p;
                    r = std::hypot(p, e_[i]);
                    e_[i + 1] = s*r;
                    s = e_[i]/r;
                    c = p/r;
                    p = c*d_[i] - s*g;
                    d_[i + 1] = h + s*(c*g + s*d_[i]);

                    // Accumulate the transformation on the first row of the
                    // rotation product alone, which is the only one the
                    // weights of the quadrature are built from
                    h = q_[i + 1];
                    q_[i + 1] = s*q_[i] + c*h;
                    q_[i] = c*q_[i] - s*h;
                }

                p = -s*s2*c3*el1*e_[l]/dl1;
                e_[l] = s*p;
                d_[l] = c*p;
            }
            while (mag(e_[l]) > eps*tst1); // Convergence check
        }

        d_[l] = d_[l] + f;
        e_[l] = 0.0;
    }

    // Sort the eigenvalues in ascending order, carrying the first row of the
    // rotation product along with them
    for (label i = 0; i < nNodes - 1; i++)
    {
        label k = i;
        scalar p = d_[i];

        for (label j = i + 1; j < nNodes; j++)
        {
            if (d_[j] < p)
            {
                k = j;
                p = d_[j];
            }
        }

        if (k != i)
        {
            d_[k] = d_[i];
            d_[i] = p;

            p = q_[i];
            q_[i] = q_[k];
            q_[k] = p;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::GolubWelsch::quadrature
(
    const UList<scalar>& alpha,
    const UList<scalar>& beta,
    const label nNodes,
    const scalar m0,
    UList<scalar>& abscissae,
    UList<scalar>& weights
)
{
    if (nNodes < 1)
    {
        FatalErrorInFunction
            << "A Gaussian quadrature needs at least one node." << nl
            << "    Number of quadrature nodes: " << nNodes << nl
            << exit(FatalError);
    }

    if (nNodes > alpha.size() || nNodes > beta.size())
    {
        FatalErrorInFunction
            << "The recurrence relationship is too short for the requested "
            << "number of quadrature nodes." << nl
            << "    Number of quadrature nodes: " << nNodes << nl
            << "    Size of the alpha list: " << alpha.size() << nl
            << "    Size of the beta list: " << beta.size() << nl
            << exit(FatalError);
    }

    if (nNodes > abscissae.size() || nNodes > weights.size())
    {
        FatalErrorInFunction
            << "The quadrature does not have room for the requested number "
            << "of nodes." << nl
            << "    Number of quadrature nodes: " << nNodes << nl
            << "    Size of the abscissae list: " << abscissae.size() << nl
            << "    Size of the weights list: " << weights.size() << nl
            << exit(FatalError);
    }

    // The coefficients whose square roots form the sub-diagonal have to be
    // positive. This is a condition on memory, and not only on the answer: a
    // NaN, which the square root of a negative coefficient produces, poisons
    // the norm the sweep measures its sub-diagonal against, so no entry of
    // the chain is ever found small, the search for a splitting point walks
    // past its end, and the sweep then reads and writes outside the work
    // arrays. The test is written as a negation so that a coefficient that
    // is already a NaN, which compares false against everything, is refused
    // here rather than carried into the sweep.
    for (label i = 1; i < nNodes; i++)
    {
        if (!(beta[i] >= 0))
        {
            FatalErrorInFunction
                << "The coefficients of the recurrence relationship of a "
                << "measure are positive, and the ones of a realizable "
                << "moment set are strictly positive." << nl
                << "    Coefficient: beta[" << i << "]" << nl
                << "    Value: " << beta[i] << nl
                << "    Number of quadrature nodes: " << nNodes << nl
                << "    Coefficients: " << beta << nl
                << exit(FatalError);
        }
    }

    if (nNodes == 1)
    {
        abscissae[0] = alpha[0];
        weights[0] = m0;

        return;
    }

    resizeWork(nNodes);

    // Jacobi matrix of the recurrence relationship, in the diagonal and
    // sub-diagonal form the QL sweep consumes. The sub-diagonal is shifted so
    // that e_[i] couples d_[i] and d_[i + 1], and its last entry closes the
    // chain.
    for (label i = 0; i < nNodes; i++)
    {
        d_[i] = alpha[i];
        q_[i] = 0.0;
    }

    for (label i = 0; i < nNodes - 1; i++)
    {
        e_[i] = Foam::sqrt(beta[i + 1]);
    }

    e_[nNodes - 1] = 0.0;

    // The Jacobi matrix is already tridiagonal, so the orthogonal matrix that
    // brings a symmetric matrix to that form is the identity here, and the
    // first row of the accumulated rotations starts as the first unit vector
    q_[0] = 1.0;

    symmTridiagonalQL(nNodes);

    // The rotations preserve the norm of q_, so the weights add up to m0
    for (label i = 0; i < nNodes; i++)
    {
        abscissae[i] = d_[i];
        weights[i] = m0*sqr(q_[i]);
    }
}


// ************************************************************************* //
