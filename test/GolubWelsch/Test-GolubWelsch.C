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

Application
    Test-GolubWelsch

Description
    Test the Gaussian quadrature built from the coefficients of the three-term
    recurrence relationship of the monic orthogonal polynomials.

    Two checks carry the weight, and neither of them records a value of a
    previous run:

    1. The classical families of orthogonal polynomials, whose nodes and
       weights are known in closed form.

    2. The moments of the Jacobi matrix. The quadrature diagonalises J, so
       \f$ \sum_i w_i x_i^k = m_0 (J^k)_{00} \f$ holds exactly, for every k and
       for any recurrence relationship with positive beta. That identity fails
       for any corruption of the eigenvalues or of the first row of the
       eigenvectors, and it covers the recurrence relationships no closed form
       is available for.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "scalarList.H"
#include "GolubWelsch.H"

#include <limits>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nChecks = 0;

//- Compare against an expected value, relative to the scale on which the
//  value is built. A moment of odd order is null by cancellation of terms
//  that are not small, so the scale of the sum, and not the value it adds up
//  to, is what says whether the answer is right.
void compareScalar
(
    const scalar computed,
    const scalar expected,
    const scalar scale,
    const scalar tolerance,
    const string& name
)
{
    nChecks++;

    if (mag(computed - expected) > tolerance*max(mag(expected), scale))
    {
        FatalErrorInFunction
            << "The value of " << name << " is not the expected one." << nl
            << "    Computed: " << computed << nl
            << "    Expected: " << expected << nl
            << "    Difference: " << mag(computed - expected) << nl
            << "    Scale: " << scale << nl
            << "    Tolerance: " << tolerance << nl
            << exit(FatalError);
    }
}


//- Moment of order momentOrder of the quadrature
Foam::scalar momentFromQuadrature
(
    const scalarList& abscissae,
    const scalarList& weights,
    const label momentOrder
)
{
    scalar moment = 0.0;

    forAll(weights, nodei)
    {
        moment += weights[nodei]*pow(abscissae[nodei], momentOrder);
    }

    return moment;
}


//- Scale of the sum the moment of order momentOrder adds up to
Foam::scalar momentScale
(
    const scalarList& abscissae,
    const scalarList& weights,
    const label momentOrder
)
{
    scalar scale = 0.0;

    forAll(weights, nodei)
    {
        scale += weights[nodei]*pow(mag(abscissae[nodei]), momentOrder);
    }

    return scale;
}


//- Moment of order momentOrder of the measure the recurrence relationship
//  belongs to, as the entry of index (0, 0) of the momentOrder-th power of the
//  Jacobi matrix, scaled by the moment of order zero
Foam::scalar momentFromRecurrence
(
    const scalarList& alpha,
    const scalarList& beta,
    const label nNodes,
    const scalar m0,
    const label momentOrder
)
{
    // Sub-diagonal of the Jacobi matrix
    scalarList e(nNodes, Zero);

    for (label i = 0; i < nNodes - 1; i++)
    {
        e[i] = Foam::sqrt(beta[i + 1]);
    }

    // Repeated products of the Jacobi matrix by the first unit vector
    scalarList v(nNodes, Zero);
    scalarList w(nNodes, Zero);

    v[0] = 1.0;

    for (label k = 0; k < momentOrder; k++)
    {
        for (label i = 0; i < nNodes; i++)
        {
            w[i] = alpha[i]*v[i];

            if (i > 0)
            {
                w[i] += e[i - 1]*v[i - 1];
            }

            if (i < nNodes - 1)
            {
                w[i] += e[i]*v[i + 1];
            }
        }

        v = w;
    }

    return m0*v[0];
}


//- Check that the quadrature reproduces the moments of its own recurrence
//  relationship, which it does for every order because it diagonalises the
//  Jacobi matrix
void testMomentsOfRecurrence
(
    const string& name,
    const scalarList& alpha,
    const scalarList& beta,
    const label nNodes,
    const scalar m0,
    const scalar tolerance
)
{
    Info<< "\nTesting the moments of the " << name
        << " recurrence relationship with " << nNodes << " nodes" << endl;

    scalarList abscissae(nNodes, Zero);
    scalarList weights(nNodes, Zero);

    GolubWelsch golubWelsch;

    golubWelsch.quadrature(alpha, beta, nNodes, m0, abscissae, weights);

    // The weights of a Gaussian quadrature are positive, and strictly so
    // when every coefficient of the recurrence relationship is: a null
    // coefficient splits the chain and leaves the measure supported on fewer
    // points than there are nodes, which the quadrature reports as a null
    // weight rather than as an error.
    bool strictlyPositive = true;

    for (label i = 1; i < nNodes; i++)
    {
        if (beta[i] <= 0.0)
        {
            strictlyPositive = false;
        }
    }

    // The abscissae are the eigenvalues of a symmetric matrix, sorted in
    // ascending order
    forAll(weights, nodei)
    {
        if
        (
            weights[nodei] < 0.0
         || (strictlyPositive && weights[nodei] <= 0.0)
        )
        {
            FatalErrorInFunction
                << "The quadrature has a weight it cannot have." << nl
                << "    Node: " << nodei << nl
                << "    Weights: " << weights << nl
                << "    Coefficients: " << beta << nl
                << exit(FatalError);
        }

        if (nodei > 0 && abscissae[nodei] < abscissae[nodei - 1])
        {
            FatalErrorInFunction
                << "The abscissae of the quadrature are not sorted." << nl
                << "    Abscissae: " << abscissae << nl
                << exit(FatalError);
        }
    }

    for (label momenti = 0; momenti < 2*nNodes; momenti++)
    {
        compareScalar
        (
            momentFromQuadrature(abscissae, weights, momenti),
            momentFromRecurrence(alpha, beta, nNodes, m0, momenti),
            momentScale(abscissae, weights, momenti),
            tolerance,
            name + " moment " + Foam::name(momenti)
        );
    }
}


//- Check that a recurrence relationship the quadrature cannot be built from
//  is refused rather than carried into the sweep. A negative coefficient
//  makes the sub-diagonal a NaN, and a NaN defeats the test that stops the
//  search for a splitting point, so the sweep would otherwise read and write
//  outside its work arrays.
void testRefused
(
    const scalarList& alpha,
    const scalarList& beta,
    const label nNodes,
    const string& name
)
{
    Info<< "\nTesting that " << name << " is refused" << endl;

    scalarList abscissae(nNodes, Zero);
    scalarList weights(nNodes, Zero);

    GolubWelsch golubWelsch;

    bool refused = false;

    const bool wasThrowing = FatalError.throwing(true);

    try
    {
        golubWelsch.quadrature
        (
            alpha, beta, nNodes, 1.0, abscissae, weights
        );
    }
    catch (const Foam::error&)
    {
        refused = true;
    }

    FatalError.throwing(wasThrowing);

    nChecks++;

    if (!refused)
    {
        FatalErrorInFunction
            << "A recurrence relationship the quadrature cannot be built "
            << "from was accepted." << nl
            << "    Case: " << name << nl
            << "    alpha: " << alpha << nl
            << "    beta: " << beta << nl
            << "    Abscissae: " << abscissae << nl
            << "    Weights: " << weights << nl
            << exit(FatalError);
    }

    Info<< "  refused, as it has to be" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<< setprecision(16);

    // A single node carries the moments of order zero and one alone, and its
    // abscissa is the first coefficient of the recurrence relationship
    {
        Info<< "\nTesting a quadrature of one node" << endl;

        const scalarList alpha({2.5});
        const scalarList beta({3.0});

        scalarList abscissae(1, Zero);
        scalarList weights(1, Zero);

        GolubWelsch golubWelsch;

        golubWelsch.quadrature(alpha, beta, 1, 7.0, abscissae, weights);

        compareScalar(abscissae[0], 2.5, 1.0, 1.0e-15, "abscissa");
        compareScalar(weights[0], 7.0, 1.0, 1.0e-15, "weight");
    }

    // Gauss-Legendre quadrature of three nodes, whose nodes and weights are
    // known in closed form. The recurrence relationship of the monic Legendre
    // polynomials on [-1, 1] has alpha_i = 0 and
    // beta_i = i^2/(4 i^2 - 1), and the measure has m0 = 2.
    {
        Info<< "\nTesting the Gauss-Legendre quadrature of three nodes"
            << endl;

        const label nNodes = 3;

        scalarList alpha(nNodes, Zero);
        scalarList beta(nNodes, Zero);

        for (label i = 1; i < nNodes; i++)
        {
            beta[i] = sqr(scalar(i))/(4.0*sqr(scalar(i)) - 1.0);
        }

        scalarList abscissae(nNodes, Zero);
        scalarList weights(nNodes, Zero);

        GolubWelsch golubWelsch;

        golubWelsch.quadrature(alpha, beta, nNodes, 2.0, abscissae, weights);

        const scalar x = Foam::sqrt(3.0/5.0);

        compareScalar(abscissae[0], -x, 1.0, 1.0e-14, "abscissa 0");
        compareScalar(abscissae[1], 0.0, 1.0, 1.0e-14, "abscissa 1");
        compareScalar(abscissae[2], x, 1.0, 1.0e-14, "abscissa 2");

        compareScalar(weights[0], 5.0/9.0, 1.0, 1.0e-14, "weight 0");
        compareScalar(weights[1], 8.0/9.0, 1.0, 1.0e-14, "weight 1");
        compareScalar(weights[2], 5.0/9.0, 1.0, 1.0e-14, "weight 2");
    }

    // Gauss-Hermite quadrature of the probabilists, alpha_i = 0, beta_i = i,
    // integrating the standard normal distribution. The moments of odd order
    // are null and the ones of even order 2k are the double factorial
    // (2k - 1)!!, which the quadrature of n nodes reproduces up to order
    // 2n - 1.
    {
        const label nNodes = 5;

        Info<< "\nTesting the moments of the standard normal distribution "
            << "with " << nNodes << " nodes" << endl;

        scalarList alpha(nNodes, Zero);
        scalarList beta(nNodes, Zero);

        for (label i = 1; i < nNodes; i++)
        {
            beta[i] = scalar(i);
        }

        scalarList abscissae(nNodes, Zero);
        scalarList weights(nNodes, Zero);

        GolubWelsch golubWelsch;

        golubWelsch.quadrature(alpha, beta, nNodes, 1.0, abscissae, weights);

        scalar doubleFactorial = 1.0;

        for (label momenti = 0; momenti < 2*nNodes; momenti++)
        {
            scalar expected = 0.0;

            if (momenti % 2 == 0)
            {
                if (momenti > 0)
                {
                    doubleFactorial *= scalar(momenti - 1);
                }

                expected = doubleFactorial;
            }

            compareScalar
            (
                momentFromQuadrature(abscissae, weights, momenti),
                expected,
                momentScale(abscissae, weights, momenti),
                1.0e-12,
                "normal moment " + Foam::name(momenti)
            );
        }
    }

    // The families of classical orthogonal polynomials, over the range of
    // node counts the quadrature methods of moments use
    for (label nNodes = 2; nNodes <= 10; nNodes++)
    {
        scalarList alpha(nNodes, Zero);
        scalarList beta(nNodes, Zero);

        for (label i = 1; i < nNodes; i++)
        {
            beta[i] = sqr(scalar(i))/(4.0*sqr(scalar(i)) - 1.0);
        }

        testMomentsOfRecurrence
        (
            "Legendre", alpha, beta, nNodes, 2.0, 1.0e-12
        );

        beta = Zero;

        for (label i = 1; i < nNodes; i++)
        {
            beta[i] = scalar(i);
        }

        testMomentsOfRecurrence
        (
            "Hermite", alpha, beta, nNodes, 1.0, 1.0e-12
        );

        for (label i = 0; i < nNodes; i++)
        {
            alpha[i] = 2.0*scalar(i) + 1.0;
        }

        for (label i = 1; i < nNodes; i++)
        {
            beta[i] = sqr(scalar(i));
        }

        testMomentsOfRecurrence
        (
            "Laguerre", alpha, beta, nNodes, 1.0, 1.0e-10
        );
    }

    // A recurrence relationship whose beta coefficients are at the scale the
    // realizability check calls degenerate. The quadrature has to stay
    // meaningful there, because that is where a population balance spends its
    // time near a Dirac delta.
    for (label exponent = 8; exponent <= 28; exponent += 4)
    {
        const label nNodes = 4;

        scalarList alpha(nNodes, Zero);
        scalarList beta(nNodes, Zero);

        for (label i = 0; i < nNodes; i++)
        {
            alpha[i] = 1.0 + 0.1*scalar(i);
        }

        for (label i = 1; i < nNodes; i++)
        {
            beta[i] = Foam::pow(10.0, -scalar(exponent));
        }

        testMomentsOfRecurrence
        (
            "vanishing beta of 1e-" + Foam::name(exponent),
            alpha,
            beta,
            nNodes,
            1.0,
            1.0e-10
        );
    }

    // A negative coefficient of the recurrence relationship, which
    // gaussLobattoMomentInversion::correctRecurrence can produce because it
    // builds the last one from a product of polynomial values whose sign is
    // not constrained. Before it was refused here it walked the sweep off the
    // end of the work arrays and returned a quadrature of NaN.
    {
        const label nNodes = 4;

        scalarList alpha(nNodes, Zero);
        scalarList beta(nNodes, Zero);

        for (label i = 0; i < nNodes; i++)
        {
            alpha[i] = 1.0 + 0.1*scalar(i);
        }

        beta[1] = 0.5;
        beta[2] = -0.25;
        beta[3] = 0.5;

        testRefused(alpha, beta, nNodes, "a negative coefficient");

        // The same on the last coefficient alone, which is the one
        // Gauss-Lobatto overwrites
        beta[2] = 0.25;
        beta[3] = -0.5;

        testRefused(alpha, beta, nNodes, "a negative last coefficient");

        // A coefficient that is already a NaN has to be refused as well: it
        // compares false against everything, so a test written the other way
        // round would let it through
        beta[3] = std::numeric_limits<scalar>::quiet_NaN();

        testRefused(alpha, beta, nNodes, "a NaN coefficient");
    }

    // A null coefficient is not an error. It splits the chain, which is what
    // the sweep does anyway once an entry is small enough, and the kernel
    // density functions of the extended quadrature method of moments
    // approach it as their parameter goes to zero.
    {
        const label nNodes = 4;

        scalarList alpha(nNodes, Zero);
        scalarList beta(nNodes, Zero);

        for (label i = 0; i < nNodes; i++)
        {
            alpha[i] = 1.0 + 0.1*scalar(i);
        }

        beta[1] = 0.5;
        beta[2] = 0.0;
        beta[3] = 0.5;

        testMomentsOfRecurrence
        (
            "a null coefficient", alpha, beta, nNodes, 1.0, 1.0e-12
        );
    }

    Info<< nl << "Checks performed: " << nChecks << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
