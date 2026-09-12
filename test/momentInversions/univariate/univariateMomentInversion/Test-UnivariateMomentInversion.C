/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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
    Test-UnivariateMomentInversion

Description
    Test univariateMomentInversion class and methods.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "IOdictionary.H"
#include "supportType.H"
#include "univariateMomentSet.H"
#include "univariateMomentInversion.H"
#include "newUnivariateMomentInversion.C"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void compareQuadrature
(
    scalarList& expectedWeights,
    scalarList& computedWeights,
    scalarList& expectedAbscissae,
    scalarList& computedAbscissae
)
{
    Info<< "\n" << endl;

    const scalar tolerance = 1.0e-14;

    Info<< "Comparing quadrature weights and abscissae with tolerance "
        << tolerance << "\n" << endl;

    if (expectedWeights.size() != computedWeights.size())
    {
        FatalErrorInFunction
            << "Weights vectors have different size: "
            << endl
            << "  Size of computed weights vector: " << computedWeights.size()
            << endl
            << "  Size of expected weights vector: " << expectedWeights.size()
            << endl
            << exit(FatalError);
    }

    if (expectedAbscissae.size() != computedAbscissae.size())
    {
        FatalErrorInFunction
            << "Abscissae vectors have different size: "
            << endl
            << "  Size of computed abscissae vector: "
            << computedAbscissae.size() << endl
            << "  Size of expected abscissae vector: "
            << expectedAbscissae.size()
            << endl
            << exit(FatalError);
    }

    forAll(expectedWeights, weightsi)
    {
        scalar magDiff
            = mag(expectedWeights[weightsi] - computedWeights[weightsi]);

        Info<< "  expectedWeights[" << weightsi << "] = "
            << expectedWeights[weightsi]
            << ", computedWeights[" << weightsi << "] = "
            << computedWeights[weightsi]
            << endl;

        if (magDiff >= tolerance)
        {
            FatalErrorInFunction
                << "Values of the quadrature weights do not match: "
                << endl
                << "  Position: " << weightsi
                << endl
                << "  Expected weight value " << expectedAbscissae[weightsi]
                << endl
                << "  Computed weight value " << computedAbscissae[weightsi]
                << endl
                << "  Magnitude of the difference " << magDiff
                << endl
                << exit(FatalError);
        }
    }

    Info<< "\nQuadrature weights match.\n" << endl;

    forAll(expectedAbscissae, abscissai)
    {
        scalar magDiff
            = mag(expectedAbscissae[abscissai] - computedAbscissae[abscissai]);

        Info<< "  expectedAbscissae[" << abscissai << "] = "
            << expectedAbscissae[abscissai]
            << ", computedAbscissae[" << abscissai << "] = "
            << computedAbscissae[abscissai]
            << endl;

        if (magDiff >= tolerance)
        {
            FatalErrorInFunction
                << "Values of the quadrature abscissae do not match: "
                << endl
                << "  Position: " << abscissai
                << endl
                << "  Expected abscissa value " << expectedAbscissae[abscissai]
                << endl
                << "  Computed abscissa value " << computedAbscissae[abscissai]
                << endl
                << "  Magnitude of the difference " << magDiff
                << endl
                << exit(FatalError);
        }
    }

    Info<< "\nQuadrature abscissae match.\n" << endl;
}

void showInputMoments(univariateMomentSet& inputMoments, string quadratureName)
{
    Info<< "\nTesting " << quadratureName << " quadrature\n" << endl;

    forAll(inputMoments, mi)
    {
        Info<< "  inputMoments[" << mi << "] = " << inputMoments[mi] << endl;
    }

    Info<< "\n" << endl;
}

void testQuadrature
(
    univariateMomentSet& m,
    dictionary& dict,
    scalarList& expectedWeights,
    scalarList& expectedAbscissae,
    string quadratureName,
    label nMaxNodes = 0
)
{
    showInputMoments(m, quadratureName);

    autoPtr<univariateMomentInversion> inversion
    (
        univariateMomentInversion::New(dict, nMaxNodes)
    );

    inversion().invert(m, 0, 1);

    scalarList weights(inversion().weights());
    scalarList abscissae(inversion().abscissae());

    Info << "Weights: " << weights << endl;
    Info << "Abscissae: " << abscissae << endl;
    Info << "Zetas: " << m.zetas() << endl;

    compareQuadrature(expectedWeights, weights, expectedAbscissae, abscissae);

    m.update(weights, abscissae);

    Info<< "\nMoments computed from quadrature\n" << endl;

    forAll(m, mi)
    {
        Info<< "  Moment " << mi << " = " << m[mi] << endl;
    }
}

// A moment vector that is a Dirac delta to within round-off, taken from the
// first time step of tutorials/pbeTransportFoam/TaylorDispersion, whose
// initial condition is the unit moment vector.
//
// The model number density functions GQMOM uses to extend the recurrence
// relationship to the additional nodes are singular for such a vector: the
// gamma coefficient divides by m0 m2 - m1^2 and the lognormal one by
// eta^(2n) - 1, both of which are null when the variance is null. The
// realizability check reports four realizable moments, because the zeta chain
// limps one step further on round-off, so the number of realizable moments
// alone does not exclude the extension. GQMOM must fall back to Gauss.
void testOddMomentCountGQMOM(dictionary& dict)
{
    // Exact moments of a gamma distribution with shape 2 and scale 1, so
    // m_n = (n + 1)!, which is a strictly realizable moment vector of odd
    // size. GQMOM has no restriction on the parity of the moment count: it
    // takes the number of regular quadrature nodes from the parity of the
    // number of realizable moments, and extends the recurrence relationship
    // in terms of that number alone.
    scalarList inputMoments({1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0});

    const label nMoments = inputMoments.size();
    const label nMaxNodes = 5;

    autoPtr<univariateMomentInversion> inversion
    (
        univariateMomentInversion::New(dict, nMaxNodes)
    );

    // The moment set has to reserve the room the quadrature asks for
    const label nAdditionalPoints =
        inversion().nAdditionalQuadraturePoints(nMoments);

    Info<< "\nTesting GQMOM with an odd number of moments\n" << endl;
    Info<< "  Number of moments = " << nMoments << endl;
    Info<< "  Number of additional quadrature points = "
        << nAdditionalPoints << endl;

    if (nAdditionalPoints != nMaxNodes - nMoments/2)
    {
        FatalErrorInFunction
            << "The number of additional quadrature points of GQMOM is not "
            << "the one an odd number of moments determines." << nl
            << "    Number of additional quadrature points: "
            << nAdditionalPoints << ", expected " << nMaxNodes - nMoments/2
            << nl << exit(FatalError);
    }

    univariateMomentSet m
    (
        inputMoments,
        supportType::RPlus,
        SMALL,
        0.0,
        nAdditionalPoints
    );

    showInputMoments(m, "GQMOM, odd number of moments");

    const label nRealizableMoments = m.nRealizableMoments(false);

    if (nRealizableMoments != nMoments)
    {
        FatalErrorInFunction
            << "This test needs a fully realizable moment vector." << nl
            << "    Number of realizable moments: " << nRealizableMoments
            << ", expected " << nMoments << nl
            << exit(FatalError);
    }

    inversion().invert(m);

    Info<< "\nVerifying the number of quadrature nodes...";

    if (inversion().nNodes() != nMaxNodes)
    {
        FatalErrorInFunction
            << "GQMOM did not extend the quadrature to the requested number "
            << "of nodes." << nl
            << "    Number of quadrature nodes: " << inversion().nNodes()
            << ", expected " << nMaxNodes << nl
            << exit(FatalError);
    }

    Info<< "OK" << endl;

    // The additional nodes leave the moments the regular quadrature
    // determines untouched. An odd number of moments determines
    // (nMoments - 1)/2 regular nodes, which preserve the moments of order
    // lower than nMoments - 1.
    const label nRegularNodes = (nMoments - 1)/2;
    const scalarList& weights(inversion().weights());
    const scalarList& abscissae(inversion().abscissae());

    const scalar tolerance = 1.0e-10;

    Info<< "\nVerifying moment conservation with tolerance " << tolerance
        << "\n" << endl;

    for (label mi = 0; mi < 2*nRegularNodes; mi++)
    {
        scalar momentFromQuadrature = 0.0;

        forAll(weights, nodei)
        {
            momentFromQuadrature += weights[nodei]*pow(abscissae[nodei], mi);
        }

        Info<< "  moment " << mi << " from quadrature = "
            << momentFromQuadrature << ", expected " << inputMoments[mi]
            << endl;

        if
        (
            mag(momentFromQuadrature - inputMoments[mi])
          > tolerance*max(mag(inputMoments[mi]), SMALL)
        )
        {
            FatalErrorInFunction
                << "The quadrature does not reproduce the moments it is "
                << "built from." << nl
                << "    Moment order: " << mi << nl
                << "    Moment from quadrature: " << momentFromQuadrature << nl
                << "    Expected: " << inputMoments[mi] << nl
                << exit(FatalError);
        }
    }

    Info<< endl;
}


void testDegenerateVarianceGQMOM
(
    dictionary& dict,
    const word& ndfType
)
{
    dict.set("ndfTypeRPlus", ndfType);

    scalarList inputMoments
    ({
        1.0000000019725894, 1.0100000019923152, 1.0201000020122384,
        1.0303015020323618, 1.040606028719355,  1.0510151359381159
    });

    // smallZeta is null by default, which is what lets the chain reach four
    // realizable moments
    univariateMomentSet m(inputMoments, supportType::RPlus, SMALL, 0.0, 1);

    showInputMoments(m, "GQMOM, " + ndfType + ", null variance");

    const label nRealizableMoments = m.nRealizableMoments(false);

    if (nRealizableMoments != 4)
    {
        FatalErrorInFunction
            << "This test needs a moment vector with four realizable "
            << "moments, otherwise the number of realizable moments alone "
            << "would exclude the extension of the recurrence." << nl
            << "    Number of realizable moments: " << nRealizableMoments
            << nl << exit(FatalError);
    }

    autoPtr<univariateMomentInversion> inversion
    (
        univariateMomentInversion::New(dict, 4)
    );

    inversion().invert(m, 0, 1);

    // Four realizable moments give two regular nodes. GQMOM must not have
    // added any.
    Info<< "\nVerifying the fall back to Gauss...";

    if (inversion().nNodes() != 2)
    {
        FatalErrorInFunction
            << "GQMOM did not fall back to Gauss on a moment vector with "
            << "null variance." << nl
            << "    Number of quadrature nodes: " << inversion().nNodes()
            << ", expected 2" << nl
            << exit(FatalError);
    }

    Info<< "OK" << endl;

    // The quadrature has to reproduce the moments it is built from
    const scalarList& weights(inversion().weights());
    const scalarList& abscissae(inversion().abscissae());

    for (label mi = 0; mi < 2; mi++)
    {
        scalar moment = 0;

        for (label nodei = 0; nodei < inversion().nNodes(); nodei++)
        {
            moment += weights[nodei]*pow(abscissae[nodei], mi);
        }

        const scalar magDiff = mag(moment - inputMoments[mi]);

        Info<< "  moment " << mi << " from quadrature = " << moment
            << ", expected " << inputMoments[mi] << endl;

        if (magDiff > 1.0e-12*mag(inputMoments[mi]))
        {
            FatalErrorInFunction
                << "The quadrature does not reproduce moment " << mi << nl
                << "    Computed: " << moment << nl
                << "    Expected: " << inputMoments[mi] << nl
                << exit(FatalError);
        }
    }

    Info<< "\n" << endl;
}


// Invert one distribution written at several scales of its abscissa, and
// require the quadrature to follow it. A population of droplets is the same
// population whether its size is written in metres, in microns or as a
// volume in cubic metres, so the weights have to come back unchanged and the
// abscissae scaled by the same factor.
//
// It holds because smallZeta, the threshold the zeta_k of the realizability
// check are compared with, defaults to zero. The zeta_k carry the units of
// the abscissa, so any positive value of it is a threshold on a size, and a
// coordinate small enough in its own units would be called degenerate
// however well spread it is. sizeCHyQMOM forced smallZeta to SMALL and had
// to normalise its size moments for exactly that reason; this test is what
// says the univariate inversion needs no such treatment.
//  The abscissae a quadrature is told to place a node at scale with the
//  coordinate as every other abscissa does, so they are given here in the
//  units of the unscaled one. Gauss-Lobatto needs both of them.
void testScaleInvariance
(
    dictionary& dict,
    const word& quadratureName,
    const label nMaxNodes = 0,
    const scalar minKnownAbscissa = 0,
    const scalar maxKnownAbscissa = 0
)
{
    // Six moments of a quadrature of three nodes, which the inversion
    // reproduces exactly
    const scalarList w({0.2, 0.5, 0.3});
    const scalarList x({0.5, 1.5, 3.0});
    const label nMoments = 6;

    Info<< "\nTesting the scale invariance of " << quadratureName.c_str()
        << "\n" << endl;

    autoPtr<univariateMomentInversion> inversion
    (
        univariateMomentInversion::New(dict, nMaxNodes)
    );

    const label nAdditionalPoints =
        inversion().nAdditionalQuadraturePoints(nMoments);

    // The quadrature of the distribution written at one scale of its
    // abscissa, copied out of the inverter so that two can be compared
    auto quadratureAt = [&](const scalar scale)
    {
        scalarList moments(nMoments, Zero);

        forAll(moments, k)
        {
            forAll(w, i)
            {
                moments[k] += w[i]*pow(scale*x[i], k);
            }
        }

        univariateMomentSet m
        (
            moments,
            supportType::RPlus,
            inversion().smallM0(),
            inversion().smallZeta(),
            nAdditionalPoints
        );

        inversion().invert
        (
            m,
            scale*minKnownAbscissa,
            scale*maxKnownAbscissa
        );

        return Tuple2<scalarList, scalarList>
        (
            inversion().weights(),
            inversion().abscissae()
        );
    };

    const Tuple2<scalarList, scalarList> reference(quadratureAt(1.0));

    Info<< "  at a scale of 1: weights " << reference.first()
        << ", abscissae " << reference.second() << endl;

    const scalar tolerance = 1.0e-12;

    // A micron written in metres, the volume of a micron-sized particle
    // written in cubic metres, and a scale as far the other way
    for (const scalar scale : {1.0e-6, 1.0e-18, 1.0e6})
    {
        const Tuple2<scalarList, scalarList> q(quadratureAt(scale));

        const scalarList& weights = q.first();
        const scalarList& abscissae = q.second();

        if
        (
            weights.size() != reference.first().size()
         || abscissae.size() != reference.second().size()
        )
        {
            FatalErrorInFunction
                << quadratureName << " builds a different number of nodes "
                << "when the abscissa is scaled by " << scale << "." << nl
                << "    Number of nodes: " << weights.size()
                << ", expected " << reference.first().size() << nl
                << exit(FatalError);
        }

        // Each list is measured against the largest of its own entries
        // rather than against each entry on its own: a quadrature is free
        // to place a node at zero, as Gauss-Lobatto does at the abscissa it
        // is told to, and such a node would otherwise be asked to come back
        // to no error at all.
        scalar weightScale = 0;
        scalar abscissaScale = 0;

        forAll(weights, nodei)
        {
            weightScale = max(weightScale, mag(reference.first()[nodei]));
            abscissaScale = max(abscissaScale, mag(reference.second()[nodei]));
        }

        scalar worstWeight = 0;
        scalar worstAbscissa = 0;

        forAll(weights, nodei)
        {
            worstWeight =
                max
                (
                    worstWeight,
                    mag(weights[nodei] - reference.first()[nodei])
                   /max(weightScale, SMALL)
                );

            worstAbscissa =
                max
                (
                    worstAbscissa,
                    mag(abscissae[nodei]/scale - reference.second()[nodei])
                   /max(abscissaScale, SMALL)
                );
        }

        Info<< "  at a scale of " << scale << ": the weights differ by "
            << worstWeight << " and the abscissae by " << worstAbscissa
            << endl;

        if (worstWeight > tolerance || worstAbscissa > tolerance)
        {
            FatalErrorInFunction
                << quadratureName << " does not invert the same distribution "
                << "to the same quadrature when the units of its abscissa "
                << "change." << nl
                << "    Scale: " << scale << nl
                << "    Weights: " << weights << nl
                << "    Weights at a scale of one: " << reference.first() << nl
                << "    Abscissae over the scale: " << abscissae/scale << nl
                << "    Abscissae at a scale of one: " << reference.second()
                << nl
                << "    Tolerance: " << tolerance << nl
                << exit(FatalError);
        }
    }

    Info<< "\n" << quadratureName.c_str()
        << " inverts the same quadrature at every scale.\n" << endl;
}


// Gauss-Lobatto fixes a node at each of two abscissae. Given an interval
// that is empty, the system that corrects the last coefficients of the
// recurrence is singular, and the inversion has to say so rather than return
// the 0/0 it used to, which surfaced later as a coefficient the Golub-Welsch
// algorithm found not to be finite.
//
// With exceptions on, any fatal error is caught, the Golub-Welsch one
// included, so the refusal is required to come from the Gauss-Lobatto
// inversion itself: otherwise the check would pass on the defect it is
// meant to catch.
void testLobattoRefusesAnEmptyInterval(dictionary& dict)
{
    Info<< "\nTesting that Gauss-Lobatto refuses an empty interval\n" << endl;

    // Six moments of a quadrature of three nodes
    const scalarList w({0.2, 0.5, 0.3});
    const scalarList x({0.5, 1.5, 3.0});
    const label nMoments = 6;

    scalarList moments(nMoments, Zero);

    forAll(moments, k)
    {
        forAll(w, i)
        {
            moments[k] += w[i]*pow(x[i], k);
        }
    }

    // Both at zero, which is what invert defaults to and what the
    // inversions that call it without abscissae pass; equal and away from
    // zero; and the wrong way round
    const List<Pair<scalar>> intervals
    ({
        Pair<scalar>(0.0, 0.0),
        Pair<scalar>(1.5, 1.5),
        Pair<scalar>(3.5, 0.0)
    });

    for (const Pair<scalar>& interval : intervals)
    {
        autoPtr<univariateMomentInversion> inversion
        (
            univariateMomentInversion::New(dict)
        );

        univariateMomentSet m
        (
            moments,
            supportType::RPlus,
            inversion().smallM0(),
            inversion().smallZeta(),
            inversion().nAdditionalQuadraturePoints(nMoments)
        );

        const bool throwing = FatalError.throwing(true);

        bool refusedHere = false;
        string refusedElsewhere;

        try
        {
            inversion().invert(m, interval.first(), interval.second());
        }
        catch (const Foam::error& err)
        {
            if (err.sourceFileName().find("gaussLobatto") != string::npos)
            {
                refusedHere = true;
            }
            else
            {
                refusedElsewhere = err.sourceFileName();
            }
        }

        FatalError.throwing(throwing);

        Info<< "  minKnownAbscissa " << interval.first()
            << ", maxKnownAbscissa " << interval.second() << ": ";

        if (refusedHere)
        {
            Info<< "refused" << endl;
            continue;
        }

        if (!refusedElsewhere.empty())
        {
            FatalErrorInFunction
                << "Gauss-Lobatto did not refuse an empty interval itself: "
                << "the inversion failed further on, in "
                << refusedElsewhere.c_str() << "." << nl
                << "    minKnownAbscissa: " << interval.first() << nl
                << "    maxKnownAbscissa: " << interval.second() << nl
                << exit(FatalError);
        }

        FatalErrorInFunction
            << "Gauss-Lobatto inverted a moment set on an empty interval "
            << "instead of refusing it." << nl
            << "    minKnownAbscissa: " << interval.first() << nl
            << "    maxKnownAbscissa: " << interval.second() << nl
            << "    Weights: " << inversion().weights() << nl
            << "    Abscissae: " << inversion().abscissae() << nl
            << exit(FatalError);
    }

    Info<< "\nGauss-Lobatto refuses an empty interval.\n" << endl;
}


int main(int argc, char *argv[])
{
    Info<< setprecision(16);

    Info<< "Testing univariateMomentInversion\n" << endl;
    Info<< "---------------------------------\n" << endl;

    Info<< "Reading quadraturePropertiesGauss\n" << endl;

    dictionary quadraturePropertiesGauss
    (
        IFstream("quadraturePropertiesGauss")()
    );

    Info<< "Reading quadraturePropertiesRadau\n" << endl;

    dictionary quadraturePropertiesRadau
    (
        IFstream("quadraturePropertiesRadau")()
    );

    Info<< "Reading quadraturePropertiesLobatto\n" << endl;

    dictionary quadraturePropertiesLobatto
    (
        IFstream("quadraturePropertiesLobatto")()
    );

    dictionary quadraturePropertiesGQMOM
    (
        IFstream("quadraturePropertiesGQMOM")()
    );

    // Test 1 - m = (1, 1, 1, 1)
    scalarList inputMoments1(4, 1.0);

    univariateMomentSet mGaussTest1
    (
        inputMoments1, supportType::RPlus, SMALL, SMALL
    );

    scalarList expectedWeightsTest1(1, 1.0);
    scalarList expectedAbscissaeTest1(1, 1.0);

    testQuadrature
    (
        mGaussTest1,
        quadraturePropertiesGauss,
        expectedWeightsTest1,
        expectedAbscissaeTest1,
        "Gauss"
    );

    // Test 2 - m_{i-1} = 1/m_i, with i = 1, ..., nMoments - 1
    scalarList inputMoments2(10);

    for (label mi = 1; mi < inputMoments2.size() + 1; mi++)
    {
        inputMoments2[mi - 1] = 1.0/scalar(mi);
    }

    univariateMomentSet mGaussTest2
    (
        inputMoments2, supportType::RPlus, SMALL, SMALL
    );

    scalarList expectedWeightsTest2(5);

    expectedWeightsTest2[0] = 0.1184634425280107;
    expectedWeightsTest2[1] = 0.2393143352499202;
    expectedWeightsTest2[2] = 0.2844444444446074;
    expectedWeightsTest2[3] = 0.2393143352495329;
    expectedWeightsTest2[4] = 0.1184634425279286;

    scalarList expectedAbscissaeTest2(5);

    expectedAbscissaeTest2[0] = 0.04691007703061419;
    expectedAbscissaeTest2[1] = 0.2307653449471704;
    expectedAbscissaeTest2[2] = 0.5000000000002638;
    expectedAbscissaeTest2[3] = 0.7692346550530993;
    expectedAbscissaeTest2[4] = 0.9530899229694018;

    testQuadrature
    (
        mGaussTest2,
        quadraturePropertiesGauss,
        expectedWeightsTest2,
        expectedAbscissaeTest2,
        "Gauss"
    );

    // Test 3 - m_{i-1} = 1/m_i, with i = 1, ..., nMoments - 1, Gauss-Lobatto
    scalarList inputMoments3(11);

    for (label mi = 1; mi < inputMoments3.size() + 1; mi++)
    {
        inputMoments3[mi - 1] = 1.0/scalar(mi);
    }

    univariateMomentSet mGaussTest3
    (
        inputMoments3, supportType::RPlus, SMALL, SMALL, 1
    );

    scalarList expectedWeightsTest3(6);

    expectedWeightsTest3[0] = 0.02777777777722189;
    expectedWeightsTest3[1] = 0.1598203766073751;
    expectedWeightsTest3[2] = 0.2426935942324191;
    expectedWeightsTest3[3] = 0.2604633915958309;
    expectedWeightsTest3[4] = 0.2084506671586522;
    expectedWeightsTest3[5] = 0.1007941926285008;

    scalarList expectedAbscissaeTest3(6);

    expectedAbscissaeTest3[0] = 0.0;
    expectedAbscissaeTest3[1] = 0.09853508579692917;
    expectedAbscissaeTest3[2] = 0.30453572664167;
    expectedAbscissaeTest3[3] = 0.562025189747392;
    expectedAbscissaeTest3[4] = 0.8019865821232577;
    expectedAbscissaeTest3[5] = 0.9601901429478156;

    testQuadrature
    (
        mGaussTest3,
        quadraturePropertiesRadau,
        expectedWeightsTest3,
        expectedAbscissaeTest3,
        "Gauss-Radau"
    );

    // Test 4 - m_{i-1} = 1/m_i, with i = 1, ..., nMoments - 1, Gauss-Radau
    scalarList inputMoments4(12);

    for (label mi = 1; mi < inputMoments4.size() + 1; mi++)
    {
        inputMoments4[mi - 1] = 1.0/scalar(mi);
    }

    univariateMomentSet mGaussTest4
    (
        inputMoments4, supportType::RPlus, SMALL, SMALL, 2
    );

    scalarList expectedWeightsTest4(7);

    expectedWeightsTest4[0] = 0.02380952380528316;
    expectedWeightsTest4[1] = 0.1384130236614941;
    expectedWeightsTest4[2] = 0.215872690592007;
    expectedWeightsTest4[3] = 0.2438095238141702;
    expectedWeightsTest4[4] = 0.2158726906203832;
    expectedWeightsTest4[5] = 0.1384130236945297;
    expectedWeightsTest4[6] = 0.02380952381213262;

    scalarList expectedAbscissaeTest4(7);

    expectedAbscissaeTest4[0] = 0.0;
    expectedAbscissaeTest4[1] = 0.08488805184722899;
    expectedAbscissaeTest4[2] = 0.2655756032333322;
    expectedAbscissaeTest4[3] = 0.499999999964492;
    expectedAbscissaeTest4[4] = 0.734424396710906;
    expectedAbscissaeTest4[5] = 0.9151119481304069;
    expectedAbscissaeTest4[6] = 1.0;

    testQuadrature
    (
        mGaussTest4,
        quadraturePropertiesLobatto,
        expectedWeightsTest4,
        expectedAbscissaeTest4,
        "Gauss-Lobatto"
    );

    // Test 5 - GQMOM on R
    scalarList inputMoments5(10);

    inputMoments5[0] = 1.0;
    inputMoments5[1] = 0.0;
    inputMoments5[2] = 1.0;
    inputMoments5[3] = 0.0;
    inputMoments5[4] = 3.0;
    inputMoments5[5] = 0.0;
    inputMoments5[6] = 15.0;
    inputMoments5[7] = 0.0;
    inputMoments5[8] = 105.0;
    inputMoments5[9] = 0.0;
    //inputMoments5[10] = 945.0;

    univariateMomentSet mGaussTest5
    (
        inputMoments5, supportType::R, SMALL, SMALL, 5
    );

    scalarList expectedWeightsTest5(10, 0);

    expectedWeightsTest5[0] = 4.310652630718267e-06;
    expectedWeightsTest5[1] = 0.0007580709343122131;
    expectedWeightsTest5[2] = 0.01911158050077029;
    expectedWeightsTest5[3] = 0.1354837029802678;
    expectedWeightsTest5[4] = 0.3446423349320191;
    expectedWeightsTest5[5] = 0.3446423349320191;
    expectedWeightsTest5[6] = 0.1354837029802678;
    expectedWeightsTest5[7] = 0.01911158050077032;
    expectedWeightsTest5[8] = 0.0007580709343122137;
    expectedWeightsTest5[9] = 4.310652630718305e-06;

    scalarList expectedAbscissaeTest5(10, 0);

    expectedAbscissaeTest5[0] = -4.859462828332314;
    expectedAbscissaeTest5[1] = -3.581823483551925;
    expectedAbscissaeTest5[2] = -2.484325841638955;
    expectedAbscissaeTest5[3] = -1.465989094391158;
    expectedAbscissaeTest5[4] = -0.4849357075154974;
    expectedAbscissaeTest5[5] = 0.4849357075154979;
    expectedAbscissaeTest5[6] = 1.465989094391158;
    expectedAbscissaeTest5[7] = 2.484325841638951;
    expectedAbscissaeTest5[8] = 3.581823483551929;
    expectedAbscissaeTest5[9] = 4.85946282833231;

    testQuadrature
    (
        mGaussTest5,
        quadraturePropertiesGQMOM,
        expectedWeightsTest5,
        expectedAbscissaeTest5,
        "GQMOM",
        10
    );


    // Test 6 - GQMOM on R+
    scalarList inputMoments6(10);

    for (label mi = 1; mi < inputMoments6.size() + 1; mi++)
    {
        inputMoments6[mi - 1] = 1.0/scalar(mi);
    }

    univariateMomentSet mGaussTest6
    (
        inputMoments6, supportType::RPlus, SMALL, SMALL, 5
    );

    scalarList expectedWeightsTest6(10, 0);

    expectedWeightsTest6[0] = 0.04676927763440855;
    expectedWeightsTest6[1] = 0.09238395831062521;
    expectedWeightsTest6[2] = 0.1338541470311804;
    expectedWeightsTest6[3] = 0.1656978021914576;
    expectedWeightsTest6[4] = 0.1735468576392944;
    expectedWeightsTest6[5] = 0.1798430425746221;
    expectedWeightsTest6[6] = 0.1249034785658432;
    expectedWeightsTest6[7] = 0.08286244866468077;
    expectedWeightsTest6[8] = 0.0001389703139075481;
    expectedWeightsTest6[9] = 1.707398051215797e-08;

    scalarList expectedAbscissaeTest6(10, 0);

    expectedAbscissaeTest6[0] = 0.0183308690766393;
    expectedAbscissaeTest6[1] = 0.09048528273499766;
    expectedAbscissaeTest6[2] = 0.2014909335053396;
    expectedAbscissaeTest6[3] = 0.3557178385194207 ;
    expectedAbscissaeTest6[4] = 0.5232653238609548;
    expectedAbscissaeTest6[5] = 0.7048225199412453;
    expectedAbscissaeTest6[6] = 0.8588367392334761;
    expectedAbscissaeTest6[7] = 0.9639042866192365;
    expectedAbscissaeTest6[8] = 1.148144194761292 ;
    expectedAbscissaeTest6[9] = 1.523890900642449;

    testQuadrature
    (
        mGaussTest6,
        quadraturePropertiesGQMOM,
        expectedWeightsTest6,
        expectedAbscissaeTest6,
        "GQMOM",
        10
    );

    // Test 7 - GQMOM on [0, 1]
    scalarList inputMoments7(10);

    inputMoments7[0] = 1.0/2.0;
    inputMoments7[1] = 3.0/10.0;
    inputMoments7[2] = 1.0/5.0;
    inputMoments7[3] = 1.0/7.0;
    inputMoments7[4] = 3.0/28.0;
    inputMoments7[5] = 1.0/12.0;
    inputMoments7[6] = 1.0/15.0;
    inputMoments7[7] = 3.0/55.0;
    inputMoments7[8] = 1.0/22.0;
    inputMoments7[9] = 1.0/26.0;

    univariateMomentSet mGaussTest7
    (
        inputMoments7, supportType::ZeroOne, SMALL, SMALL, 5
    );

    scalarList expectedWeightsTest7(10, 0);

    expectedWeightsTest7[0] = 0.001943873229307942;
    expectedWeightsTest7[1] = 0.007902072129949747;
    expectedWeightsTest7[2] = 0.02732546156051659;
    expectedWeightsTest7[3] = 0.05761179944540628;
    expectedWeightsTest7[4] = 0.09028841825500214;
    expectedWeightsTest7[5] = 0.10986807966067;
    expectedWeightsTest7[6] = 0.09845777421751999;
    expectedWeightsTest7[7] = 0.06927178902222311;
    expectedWeightsTest7[8] = 0.02951883479525387;
    expectedWeightsTest7[9] = 0.007811897684150718;

    scalarList expectedAbscissaeTest7(10, 0);

    expectedAbscissaeTest7[0] = 0.06964594081013871;
    expectedAbscissaeTest7[1] = 0.14017230848925;
    expectedAbscissaeTest7[2] = 0.2355328726617618;
    expectedAbscissaeTest7[3] = 0.3490332042621505;
    expectedAbscissaeTest7[4] = 0.4750293831698755;
    expectedAbscissaeTest7[5] = 0.6015400657359764;
    expectedAbscissaeTest7[6] = 0.723885477138646;
    expectedAbscissaeTest7[7] = 0.8284968176504827;
    expectedAbscissaeTest7[8] = 0.9147863479347316;
    expectedAbscissaeTest7[9] = 0.9684906429192973;

    testQuadrature
    (
        mGaussTest7,
        quadraturePropertiesGQMOM,
        expectedWeightsTest7,
        expectedAbscissaeTest7,
        "GQMOM",
        10
    );

    scalarList inputMoments8(6);

    inputMoments8[0] = 1.0;
    inputMoments8[1] = 1.13;
    inputMoments8[2] = 1.294;
    inputMoments8[3] = 1.5;
    inputMoments8[4] = 1.760;
    inputMoments8[5] = 2.090237;

    univariateMomentSet mGaussTest8
    (
        inputMoments8, supportType::RPlus, SMALL, SMALL, 7
    );

    scalarList expectedWeightsTest8(10, 0);

    expectedWeightsTest8[0] = 0.001528021029486816;
    expectedWeightsTest8[1] = 0.001386763449248363;
    expectedWeightsTest8[2] = 0.001866128388630542;
    expectedWeightsTest8[3] = 0.003630846785146271;
    expectedWeightsTest8[4] = 0.01421509421173688;
    expectedWeightsTest8[5] = 0.6853246975604294;
    expectedWeightsTest8[6] = 0.02431175018572058;
    expectedWeightsTest8[7] = 0.2675675885682817;
    expectedWeightsTest8[8] = 0.0001691025284614494;
    expectedWeightsTest8[9] = 7.292857634940761e-09;

    scalarList expectedAbscissaeTest8(10, 0);

    expectedAbscissaeTest8[0] = 0.1840361378198212;
    expectedAbscissaeTest8[1] = 0.3716851294787171;
    expectedAbscissaeTest8[2] = 0.541482301270592;
    expectedAbscissaeTest8[3] = 0.7215189555042092;
    expectedAbscissaeTest8[4] = 0.9225265146913547;
    expectedAbscissaeTest8[5] = 1.064838394996712;
    expectedAbscissaeTest8[6] = 1.153156632425781;
    expectedAbscissaeTest8[7] = 1.324616410466969;
    expectedAbscissaeTest8[8] = 1.41694749651431;
    expectedAbscissaeTest8[9] = 1.763388388042151;

    testQuadrature
    (
        mGaussTest8,
        quadraturePropertiesGQMOM,
        expectedWeightsTest8,
        expectedAbscissaeTest8,
        "GQMOM",
        10
    );

    testOddMomentCountGQMOM(quadraturePropertiesGQMOM);

    testDegenerateVarianceGQMOM(quadraturePropertiesGQMOM, "gamma");
    testDegenerateVarianceGQMOM(quadraturePropertiesGQMOM, "lognormal");

    testScaleInvariance(quadraturePropertiesGauss, "Gauss");
    testScaleInvariance(quadraturePropertiesRadau, "Gauss-Radau");
    testScaleInvariance(quadraturePropertiesLobatto, "Gauss-Lobatto", 0, 0, 3.5);

    testLobattoRefusesAnEmptyInterval(quadraturePropertiesLobatto);
    testScaleInvariance(quadraturePropertiesGQMOM, "GQMOM", 5);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
