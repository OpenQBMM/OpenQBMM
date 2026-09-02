/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2026 Alberto Passalacqua
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
    Test-UnivariateMomentSet

Description
    Test univariateMomentSet class and methods.

    Each test verifies the complete state produced by the realizability check:
    the zeta_k values, the canonical moments (support [0, 1] only), the alpha
    and beta coefficients of the recurrence relationship, and every
    realizability flag. The sizes of the alpha and beta lists are part of the
    expected values, because they are what guarantees that the recurrence
    relationship is built without writing past the end of the lists.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "IOdictionary.H"
#include "supportType.H"
#include "univariateMomentSet.H"
#include "univariateMomentInversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Tolerance used to compare computed and expected values.
//  The expected values are the ones the algorithms produce, written with the
//  17 significant digits a double needs to be recovered exactly, so every
//  comparison of this test is satisfied with a null difference. The tolerance
//  is only there to absorb the round-off a different compiler, or a different
//  optimization level, may introduce.
static const scalar testTol = SMALL;

//- Compare two lists of scalars, aborting if sizes or values differ
void compareScalarLists
(
    const string& name,
    const scalarList& computed,
    const scalarList& expected
)
{
    Info<< "\nVerifying " << name << nl << endl;

    if (computed.size() != expected.size())
    {
        FatalErrorInFunction
            << "Lists of " << name << " have different size:" << nl
            << "    Size of computed list: " << computed.size() << nl
            << "    Size of expected list: " << expected.size() << nl
            << "    Computed: " << computed << nl
            << "    Expected: " << expected << nl
            << exit(FatalError);
    }

    forAll(computed, i)
    {
        Info<< "  expected[" << i << "] = " << setprecision(17)
            << expected[i] << ", computed[" << i << "] = "
            << computed[i] << endl;
    }

    forAll(computed, i)
    {
        const scalar magDiff = mag(computed[i] - expected[i]);

        if (magDiff >= testTol)
        {
            FatalErrorInFunction
                << "Values of " << name << " do not match:" << nl
                << "    Position: " << i << nl
                << "    Expected value: " << setprecision(17) << expected[i]
                << nl
                << "    Computed value: " << computed[i] << nl
                << "    Magnitude of the difference: " << magDiff << nl
                << "    Tolerance: " << testTol << nl
                << exit(FatalError);
        }
    }

    Info<< "\nValues of " << name << " match.\n" << endl;
}


//- Compare a label, aborting if the values differ
void compareLabel
(
    const string& name,
    const label computed,
    const label expected
)
{
    Info<< "Verifying " << name << " (expected " << expected
        << ", computed " << computed << ")...";

    if (computed != expected)
    {
        FatalErrorInFunction
            << "Value of " << name << " does not match:" << nl
            << "    Expected value: " << expected << nl
            << "    Computed value: " << computed << nl
            << exit(FatalError);
    }

    Info<< "OK" << endl;
}


//- Compare a bool, aborting if the values differ
void compareFlag
(
    const string& name,
    const bool computed,
    const bool expected
)
{
    Info<< "Verifying " << name << " (expected " << expected
        << ", computed " << computed << ")...";

    if (computed != expected)
    {
        FatalErrorInFunction
            << "Value of " << name << " does not match:" << nl
            << "    Expected value: " << expected << nl
            << "    Computed value: " << computed << nl
            << exit(FatalError);
    }

    Info<< "OK" << endl;
}


//- Verify the full realizability status of a moment set
void compareStatus
(
    univariateMomentSet& moments,
    const label nRealizableMoments,
    const bool fullyRealizable,
    const bool subsetRealizable,
    const bool onMomentSpaceBoundary,
    const bool degenerate
)
{
    Info<< "\nVerifying realizability status\n" << endl;

    compareLabel
    (
        "nRealizableMoments", moments.nRealizableMoments(false),
        nRealizableMoments
    );

    compareFlag
    (
        "isFullyRealizable", moments.isFullyRealizable(false), fullyRealizable
    );

    compareFlag
    (
        "isSubsetRealizable", moments.isSubsetRealizable(false),
        subsetRealizable
    );

    compareFlag
    (
        "isOnMomentSpaceBoundary", moments.isOnMomentSpaceBoundary(false),
        onMomentSpaceBoundary
    );

    compareFlag("isDegenerate", moments.isDegenerate(), degenerate);
}


//- Report the moment vector under test
void showInputMoments(const string& name, const scalarList& m)
{
    Info<< nl << "Testing " << name << nl
        << "----------------------------------------" << endl;

    forAll(m, mi)
    {
        Info<< "  inputMoments[" << mi << "] = " << setprecision(17)
            << m[mi] << endl;
    }
}


//- Build the moment vector of a Gaussian distribution with the given mean
//  and standard deviation
scalarList gaussianMoments
(
    const label nMoments,
    const scalar mu,
    const scalar sigma
)
{
    scalarList m(nMoments, Zero);

    // Central moments of a Gaussian: mu_{2k} = (2k - 1)!! sigma^{2k},
    // mu_{2k+1} = 0. Raw moments follow from the binomial expansion.
    scalarList centralMoments(nMoments, Zero);
    centralMoments[0] = 1.0;

    for (label k = 2; k < nMoments; k++)
    {
        if (k % 2 == 0)
        {
            centralMoments[k] = scalar(k - 1)*sqr(sigma)*centralMoments[k - 2];
        }
    }

    for (label k = 0; k < nMoments; k++)
    {
        scalar binomial = 1.0;

        for (label j = 0; j <= k; j++)
        {
            m[k] += binomial*centralMoments[j]*pow(mu, k - j);

            binomial *= scalar(k - j)/scalar(j + 1);
        }
    }

    return m;
}


//- Build the moment vector m_k = 1/(k + 1), i.e. the moments of a uniform
//  distribution over [0, 1]
scalarList uniformMoments(const label nMoments)
{
    scalarList m(nMoments);

    forAll(m, k)
    {
        m[k] = 1.0/scalar(k + 1);
    }

    return m;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Moments of a Gaussian with mu = -2 and sigma = 1, over R.
//
// The recurrence relationship of the monic orthogonal polynomials of a
// Gaussian measure is known analytically - alpha_k = mu and beta_k = k sigma^2,
// with beta_0 = m0 - so the expected values are exact rather than recorded
// from a previous run.
//
// The odd number of moments is what this test is for: the walk of the
// recurrence has to stop at the last coefficient the moments determine.
// Continuing it by one would write alpha_[nD], with nD = (nMoments - 1)/2,
// which is both past the end of the alpha list and not determined by the
// moments.
void testGaussianMomentsR()
{
    const label nMoments = 9;
    const scalar mu = -2;
    const scalar sigma = 1;

    scalarList inputMoments(gaussianMoments(nMoments, mu, sigma));

    showInputMoments
    (
        "Gaussian moments, support R, odd moment count", inputMoments
    );

    univariateMomentSet moments(inputMoments, supportType::R, SMALL, SMALL);

    compareStatus(moments, nMoments, true, true, false, false);

    // With nMoments = 2p + 1 the moments determine alpha_0, ..., alpha_{p-1}
    // and beta_0, ..., beta_p: alpha_p would require the moment of order
    // 2p + 1, which is not available.
    scalarList expectedAlpha(4, mu);

    // beta_0 = m0, beta_k = k sigma^2
    scalarList expectedBeta(5);
    expectedBeta[0] = inputMoments[0];
    for (label k = 1; k < expectedBeta.size(); k++)
    {
        expectedBeta[k] = scalar(k)*sqr(sigma);
    }

    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// The same Gaussian, with an even number of moments. This is the control for
// the case above: the last alpha written by the recurrence is alpha_[nD] with
// nD = nMoments/2 - 1, one position lower than for an odd moment count.
void testGaussianMomentsREvenCount()
{
    const label nMoments = 8;
    const scalar mu = -2;
    const scalar sigma = 1;

    scalarList inputMoments(gaussianMoments(nMoments, mu, sigma));

    showInputMoments
    (
        "Gaussian moments, support R, even moment count", inputMoments
    );

    univariateMomentSet moments(inputMoments, supportType::R, SMALL, SMALL);

    compareStatus(moments, nMoments, true, true, false, false);

    scalarList expectedAlpha(4, mu);

    // The last entry of beta is not used with an even number of moments
    scalarList expectedBeta({1.0, 1.0, 2.0, 3.0, 0.0});

    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// A fully realizable moment vector over R+, built as m_k = 1/(k + 1). These
// are the moments of a uniform distribution over [0, 1], whose recurrence
// relationship is alpha_k = 1/2, beta_k = k^2/(4(4k^2 - 1)). With 20 moments
// the higher coefficients lose accuracy, so the expected values are the ones
// the algorithm produces rather than the analytical ones.
void testFullyRealizableMomentVectorRPlus()
{
    const label nMoments = 20;

    scalarList inputMoments(uniformMoments(nMoments));

    showInputMoments
    (
        "fully realizable moment vector, support R+", inputMoments
    );

    univariateMomentSet moments(inputMoments, supportType::RPlus, SMALL, SMALL);

    compareStatus(moments, nMoments, true, true, false, false);

    scalarList expectedZetas
    ({
        0.5,
        0.16666666666666663,
        0.3333333333333338,
        0.19999999999999904,
        0.2999999999999967,
        0.21428571428575172,
        0.2857142857140651,
        0.22222222222278568,
        0.2777777777779511,
        0.22727272725651565,
        0.2727272728794396,
        0.23076922999164143,
        0.26923077241136395,
        0.23333332950536528,
        0.26666660906689227,
        0.23529475185140633,
        0.26470105237163133,
        0.23686837764181434,
        0.2630193512153205
    });

    scalarList expectedAlpha
    ({
        0.5,
        0.5000000000000004,
        0.4999999999999958,
        0.4999999999998168,
        0.5000000000007367,
        0.5000000001359552,
        0.5000000024030054,
        0.49999993857225755,
        0.4999958042230377,
        0.49988772885713484
    });

    scalarList expectedBeta
    ({
        1.0,
        0.08333333333333331,
        0.06666666666666644,
        0.06428571428572481,
        0.06349206349217545,
        0.06313131312684929,
        0.06293706276010855,
        0.06282051253204479,
        0.0627452536074504,
        0.06269930883534924,
        0.0
    });

    compareScalarLists("zetas", moments.zetas(), expectedZetas);
    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// The same moment vector as above, with m_3 set to zero. Only the first three
// moments are realizable, and the recurrence stops at zeta_2.
void testSubsetRealizableMomentVectorRPlus()
{
    const label nMoments = 10;

    scalarList inputMoments(uniformMoments(nMoments));
    inputMoments[3] = 0.0;

    showInputMoments
    (
        "moment vector with three realizable moments, support R+", inputMoments
    );

    univariateMomentSet moments(inputMoments, supportType::RPlus, SMALL, SMALL);

    compareStatus(moments, 3, false, true, false, false);

    scalarList expectedZetas(nMoments - 1, Zero);
    expectedZetas[0] = 0.5;
    expectedZetas[1] = 0.16666666666666663;
    expectedZetas[2] = -2.666666666666667;

    scalarList expectedAlpha(5, Zero);
    expectedAlpha[0] = 0.5;
    expectedAlpha[1] = -2.5000000000000004;

    scalarList expectedBeta(6, Zero);
    expectedBeta[0] = 1.0;
    expectedBeta[1] = 0.08333333333333331;

    compareScalarLists("zetas", moments.zetas(), expectedZetas);
    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// The unit moment vector m = (1 1 1 ... 1) is the moment vector of a Dirac
// delta in x = 1. Only two moments are realizable.
void testUnitMomentVectorRPlus()
{
    const label nMoments = 10;

    scalarList inputMoments(nMoments, scalar(1));

    showInputMoments("unit moment vector, support R+", inputMoments);

    univariateMomentSet moments(inputMoments, supportType::RPlus, SMALL, SMALL);

    compareStatus(moments, 2, false, true, false, false);

    scalarList expectedZetas(nMoments - 1, Zero);
    expectedZetas[0] = 1.0;

    scalarList expectedAlpha(5, Zero);
    expectedAlpha[0] = 1.0;

    scalarList expectedBeta(6, Zero);
    expectedBeta[0] = 1.0;

    compareScalarLists("zetas", moments.zetas(), expectedZetas);
    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// The moments m_k = 1/(k + 1) with support [0, 1]. All the canonical moments
// belong to [0, 1], so the moment vector is fully realizable.
void testRealizableCanonicalMoments()
{
    const label nMoments = 9;

    scalarList inputMoments(uniformMoments(nMoments));

    showInputMoments
    (
        "realizable canonical moments, support [0, 1]", inputMoments
    );

    univariateMomentSet moments
    (
        inputMoments, supportType::ZeroOne, SMALL, SMALL
    );

    compareStatus(moments, nMoments, true, true, false, false);

    scalarList expectedZetas
    ({
        0.5,
        0.16666666666666663,
        0.3333333333333338,
        0.19999999999999904,
        0.2999999999999967,
        0.21428571428575172,
        0.2857142857140651,
        0.22222222222278568
    });

    scalarList expectedCanonicalMoments
    ({
        0.5,
        0.33333333333333326,
        0.5000000000000007,
        0.39999999999999863,
        0.49999999999999334,
        0.4285714285714977,
        0.49999999999967437,
        0.4444444444452819
    });

    scalarList expectedAlpha
    ({
        0.5,
        0.5000000000000004,
        0.4999999999999958,
        0.4999999999998168
    });

    scalarList expectedBeta
    ({
        1.0,
        0.08333333333333331,
        0.06666666666666644,
        0.06428571428572481,
        0.06349206349217545
    });

    compareScalarLists("zetas", moments.zetas(), expectedZetas);
    compareScalarLists
    (
        "canonical moments", moments.canonicalMoments(),
        expectedCanonicalMoments
    );
    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// The same moments with support [0, 1] and an even number of moments, which
// is the usual configuration of QMOM. The moment vector is fully realizable,
// and isFullyRealizable() has to agree with nRealizableMoments().
void testRealizableCanonicalMomentsEvenCount()
{
    const label nMoments = 10;

    scalarList inputMoments(uniformMoments(nMoments));

    showInputMoments
    (
        "realizable canonical moments, support [0, 1], even moment count",
        inputMoments
    );

    univariateMomentSet moments
    (
        inputMoments, supportType::ZeroOne, SMALL, SMALL
    );

    compareStatus(moments, nMoments, true, true, false, false);

    scalarList expectedCanonicalMoments
    ({
        0.5,
        0.33333333333333326,
        0.5000000000000007,
        0.39999999999999863,
        0.49999999999999334,
        0.4285714285714977,
        0.49999999999967437,
        0.4444444444452819,
        0.5000000000010657
    });

    scalarList expectedAlpha
    ({
        0.5,
        0.5000000000000004,
        0.4999999999999958,
        0.4999999999998168,
        0.5000000000007367
    });

    scalarList expectedBeta
    ({
        1.0,
        0.08333333333333331,
        0.06666666666666644,
        0.06428571428572481,
        0.06349206349217545,
        0.0
    });

    compareScalarLists
    (
        "canonical moments", moments.canonicalMoments(),
        expectedCanonicalMoments
    );
    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// The same moment vector with m_2 = -1. The second canonical moment leaves
// [0, 1], so only two moments are realizable.
void testUnrealizableCanonicalMoments()
{
    const label nMoments = 9;

    scalarList inputMoments(uniformMoments(nMoments));
    inputMoments[2] = -1.0;

    showInputMoments
    (
        "unrealizable canonical moments, support [0, 1]", inputMoments
    );

    univariateMomentSet moments
    (
        inputMoments, supportType::ZeroOne, SMALL, SMALL
    );

    compareStatus(moments, 2, false, true, false, false);

    scalarList expectedZetas(nMoments - 1, Zero);
    expectedZetas[0] = 0.5;
    expectedZetas[1] = -2.5;

    scalarList expectedCanonicalMoments(nMoments - 1, Zero);
    expectedCanonicalMoments[0] = 0.5;
    expectedCanonicalMoments[1] = -5.0;

    scalarList expectedAlpha(4, Zero);
    expectedAlpha[0] = 0.5;

    scalarList expectedBeta(5, Zero);
    expectedBeta[0] = 1.0;
    expectedBeta[1] = -1.25;

    compareScalarLists("zetas", moments.zetas(), expectedZetas);
    compareScalarLists
    (
        "canonical moments", moments.canonicalMoments(),
        expectedCanonicalMoments
    );
    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// Additional quadrature points reserve room in the recurrence relationship,
// and in the zeta and canonical moment lists, for the extra nodes added by
// Gauss-Radau, Gauss-Lobatto and GQMOM. The additional entries are not
// touched by the realizability check and must be left at zero.
void testAdditionalQuadraturePoints()
{
    const label nMoments = 10;
    const label nAdditionalQuadraturePoints = 2;

    scalarList inputMoments(uniformMoments(nMoments));

    showInputMoments
    (
        "moment vector with additional quadrature points, support R+",
        inputMoments
    );

    univariateMomentSet moments
    (
        inputMoments,
        supportType::RPlus,
        SMALL,
        SMALL,
        nAdditionalQuadraturePoints
    );

    compareStatus(moments, nMoments, true, true, false, false);

    // Room for 2*(nMoments/2 + nAdditionalQuadraturePoints) - 1 zetas
    scalarList expectedZetas(13, Zero);
    expectedZetas[0] = 0.5;
    expectedZetas[1] = 0.16666666666666663;
    expectedZetas[2] = 0.3333333333333338;
    expectedZetas[3] = 0.19999999999999904;
    expectedZetas[4] = 0.2999999999999967;
    expectedZetas[5] = 0.21428571428575172;
    expectedZetas[6] = 0.2857142857140651;
    expectedZetas[7] = 0.22222222222278568;
    expectedZetas[8] = 0.2777777777779511;

    scalarList expectedAlpha(7, Zero);
    expectedAlpha[0] = 0.5;
    expectedAlpha[1] = 0.5000000000000004;
    expectedAlpha[2] = 0.4999999999999958;
    expectedAlpha[3] = 0.4999999999998168;
    expectedAlpha[4] = 0.5000000000007367;

    scalarList expectedBeta(8, Zero);
    expectedBeta[0] = 1.0;
    expectedBeta[1] = 0.08333333333333331;
    expectedBeta[2] = 0.06666666666666644;
    expectedBeta[3] = 0.06428571428572481;
    expectedBeta[4] = 0.06349206349217545;

    compareScalarLists("zetas", moments.zetas(), expectedZetas);
    compareScalarLists("alpha", moments.alphaRecurrence(), expectedAlpha);
    compareScalarLists("beta", moments.betaRecurrence(), expectedBeta);
}


// A moment vector with only two moments is realizable for every support, as
// long as the first-order moment is compatible with the support. The
// recurrence relationship of the single-node quadrature must still be set.
void testTwoMomentVectors()
{
    scalarList inputMoments({scalar(1), scalar(0.25)});

    // Support [0, 1]
    {
        showInputMoments("two moments, support [0, 1]", inputMoments);

        univariateMomentSet moments
        (
            inputMoments, supportType::ZeroOne, SMALL, SMALL
        );

        compareStatus(moments, 2, true, true, false, false);

        compareScalarLists
        (
            "zetas", moments.zetas(), scalarList(1, scalar(0.25))
        );
        compareScalarLists
        (
            "canonical moments", moments.canonicalMoments(),
            scalarList(1, scalar(0.25))
        );
        compareScalarLists
        (
            "alpha", moments.alphaRecurrence(), scalarList(1, scalar(0.25))
        );
        compareScalarLists
        (
            "beta", moments.betaRecurrence(), scalarList({1.0, 0.0})
        );
    }

    // Support R+
    {
        showInputMoments("two moments, support R+", inputMoments);

        univariateMomentSet moments
        (
            inputMoments, supportType::RPlus, SMALL, SMALL
        );

        compareStatus(moments, 2, true, true, false, false);

        compareScalarLists
        (
            "zetas", moments.zetas(), scalarList(1, scalar(0.25))
        );
        compareScalarLists
        (
            "alpha", moments.alphaRecurrence(), scalarList(1, scalar(0.25))
        );
        compareScalarLists
        (
            "beta", moments.betaRecurrence(), scalarList({1.0, 0.0})
        );
    }

    // Support R, with a negative first-order moment
    {
        scalarList negativeMeanMoments({scalar(1), scalar(-3)});

        showInputMoments("two moments, support R", negativeMeanMoments);

        univariateMomentSet moments
        (
            negativeMeanMoments, supportType::R, SMALL, SMALL
        );

        compareStatus(moments, 2, true, true, false, false);

        compareScalarLists("zetas", moments.zetas(), scalarList(1, scalar(-3)));
        compareScalarLists
        (
            "alpha", moments.alphaRecurrence(), scalarList(1, scalar(-3))
        );
        compareScalarLists
        (
            "beta", moments.betaRecurrence(), scalarList({1.0, 0.0})
        );
    }
}


// For measures with support over R+ and [0, 1] the zeta chain is the primary
// representation of the recurrence relationship, and zetasToRecurrence()
// recovers alpha and beta from it. Recovering them from the chain the
// realizability check just produced has to return the coefficients the check
// computed, which is a round trip through
//
//     zeta_{2j-1} = beta_j/zeta_{2j-2},  zeta_{2j} = alpha_j - zeta_{2j-1}
//
// An odd number of moments leaves an even number of zeta_k, the case in which
// one more beta than alpha is determined by the chain: beta_q needs only
// zeta_{2q-1}, while alpha_q would need zeta_{2q}, which is not available. The
// last beta must still be recovered, because a null beta decouples the last
// row of the Jacobi matrix and adds a spurious node of null weight.
void testZetasToRecurrence()
{
    const label nMoments = 9;

    scalarList inputMoments(uniformMoments(nMoments));

    showInputMoments
    (
        "recurrence relationship recovered from the zeta_k, support R+",
        inputMoments
    );

    univariateMomentSet moments(inputMoments, supportType::RPlus, SMALL, SMALL);

    compareStatus(moments, nMoments, true, true, false, false);

    // The coefficients built by the realizability check
    scalarList expectedAlpha(moments.alphaRecurrence());
    scalarList expectedBeta(moments.betaRecurrence());

    // An odd number of moments leaves an even number of zeta_k
    const label nZeta = nMoments - 1;

    if (nZeta % 2 != 0)
    {
        FatalErrorInFunction
            << "This test requires an even number of zeta_k." << nl
            << "    Number of zeta_k: " << nZeta << nl
            << exit(FatalError);
    }

    // The last beta must be non-null, otherwise the test cannot tell a
    // recovered coefficient from a missing one
    if (mag(expectedBeta[nZeta/2]) <= SMALL)
    {
        FatalErrorInFunction
            << "The last beta coefficient is null, so the test cannot detect "
            << "a coefficient that is not recovered." << nl
            << "    beta: " << expectedBeta << nl
            << exit(FatalError);
    }

    // The coefficients are recovered into lists of the caller, which start
    // null so that a coefficient zetasToRecurrence does not write is visible
    scalarList recoveredAlpha(expectedAlpha.size(), Zero);
    scalarList recoveredBeta(expectedBeta.size(), Zero);

    moments.zetasToRecurrence(nZeta, recoveredAlpha, recoveredBeta);

    compareScalarLists("alpha", recoveredAlpha, expectedAlpha);
    compareScalarLists("beta", recoveredBeta, expectedBeta);

    // The coefficients of the moment set must not have been altered
    compareScalarLists
    (
        "alpha of the moment set", moments.alphaRecurrence(), expectedAlpha
    );
    compareScalarLists
    (
        "beta of the moment set", moments.betaRecurrence(), expectedBeta
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<< "Testing univariateMomentSet\n" << endl;

    testGaussianMomentsR();
    testGaussianMomentsREvenCount();
    testFullyRealizableMomentVectorRPlus();
    testSubsetRealizableMomentVectorRPlus();
    testUnitMomentVectorRPlus();
    testRealizableCanonicalMoments();
    testRealizableCanonicalMomentsEvenCount();
    testUnrealizableCanonicalMoments();
    testAdditionalQuadraturePoints();
    testTwoMomentVectors();
    testZetasToRecurrence();

    Info<< "\nAll tests passed.\n" << endl;
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
