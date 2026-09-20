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
    scalarList& computedAbscissae,
    const scalar tolerance = 1.0e-14
)
{
    Info<< "\n" << endl;

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
    label nMaxNodes = 0,
    const scalar tolerance = 1.0e-14
)
{
    showInputMoments(m, quadratureName);

    autoPtr<univariateMomentInversion> inversion
    (
        univariateMomentInversion::New(dict, nMaxNodes)
    );

    inversion().invert(m);

    scalarList weights(inversion().weights());
    scalarList abscissae(inversion().abscissae());

    Info << "Weights: " << weights << endl;
    Info << "Abscissae: " << abscissae << endl;
    Info << "Zetas: " << m.zetas() << endl;

    compareQuadrature
    (
        expectedWeights,
        weights,
        expectedAbscissae,
        abscissae,
        tolerance
    );

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

    // The additional nodes leave every moment the quadrature is given
    // untouched. An odd number of moments determines (nMoments - 1)/2
    // regular nodes, which preserve the moments of order lower than
    // nMoments - 1, and the last moment determines the coefficient the
    // closure of the recurrence relationship starts from, so that it is
    // preserved too.
    const scalarList& weights(inversion().weights());
    const scalarList& abscissae(inversion().abscissae());

    const scalar tolerance = 1.0e-10;

    Info<< "\nVerifying moment conservation with tolerance " << tolerance
        << "\n" << endl;

    for (label mi = 0; mi < nMoments; mi++)
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


//- The generalized quadrature continues the recurrence relationship of the
//  moments with the one of the polynomials of a model distribution, so that
//  given the moments of that distribution it has to return its quadrature:
//  the moments it was given, and the ones of higher order the additional
//  nodes reach, up to 2 nNodes - 1. This holds whether it is given 2n moments,
//  as a case with the moments of n regular nodes does, or 2n + 1, which is
//  the case the method is written for, and for any number of additional
//  nodes.
void testModelDistributionsGQMOM()
{
    Info<< "\nTesting that GQMOM reproduces its model distributions\n"
        << endl;

    struct Model
    {
        word name;
        supportType support;
        string entries;
        scalarList moments;
    };

    const label nExact = 14;

    List<Model> models(4);

    // Gaussian of mean 0.3 and standard deviation 1.2
    models[0] = {"Gaussian", supportType::R, "nu 1;", scalarList(nExact)};
    models[0].moments[0] = 1;
    models[0].moments[1] = 0.3;

    // Gamma of shape 3 and scale 0.5, lognormal of mu = 0 and sigma = 0.3,
    // beta(2, 5)
    models[1] =
        {"gamma", supportType::RPlus, "ndfTypeRPlus gamma;", scalarList(nExact)};
    models[2] =
    {
        "lognormal", supportType::RPlus, "ndfTypeRPlus lognormal;",
        scalarList(nExact)
    };
    models[3] = {"beta", supportType::ZeroOne, "", scalarList(nExact)};

    models[1].moments[0] = 1;
    models[3].moments[0] = 1;

    for (label k = 0; k < nExact; k++)
    {
        if (k > 1)
        {
            models[0].moments[k] =
                0.3*models[0].moments[k - 1]
              + (k - 1)*sqr(1.2)*models[0].moments[k - 2];
        }

        if (k > 0)
        {
            models[1].moments[k] = models[1].moments[k - 1]*(3.0 + k - 1)*0.5;

            models[3].moments[k] =
                models[3].moments[k - 1]*(2.0 + k - 1)/(7.0 + k - 1);
        }

        models[2].moments[k] = Foam::exp(0.5*sqr(0.3*k));
    }

    for (const Model& model : models)
    {
        for (const label nMoments : {6, 7})
        {
            for (const label nAdditional : {1, 3})
            {
                const label nNodes = nMoments/2 + nAdditional;

                dictionary dict
                (
                    IStringStream
                    (
                        "univariateMomentInversion GQMOM; " + model.entries
                    )()
                );

                autoPtr<univariateMomentInversion> inversion
                (
                    univariateMomentInversion::New(dict, nNodes)
                );

                scalarList inputMoments(nMoments);

                forAll(inputMoments, mi)
                {
                    inputMoments[mi] = model.moments[mi];
                }

                univariateMomentSet m
                (
                    inputMoments,
                    model.support,
                    SMALL,
                    SMALL,
                    inversion().nAdditionalQuadraturePoints(nMoments)
                );

                inversion().invert(m);

                const scalarList& weights(inversion().weights());
                const scalarList& abscissae(inversion().abscissae());

                if (inversion().nNodes() != nNodes)
                {
                    FatalErrorInFunction
                        << "GQMOM did not extend the quadrature of the "
                        << model.name << " distribution." << nl
                        << "    Number of quadrature nodes: "
                        << inversion().nNodes() << ", expected " << nNodes
                        << nl << exit(FatalError);
                }

                scalar largest = 0;

                for (label mi = 0; mi < 2*nNodes; mi++)
                {
                    scalar momentFromQuadrature = 0.0;

                    forAll(weights, nodei)
                    {
                        momentFromQuadrature +=
                            weights[nodei]*pow(abscissae[nodei], mi);
                    }

                    // The first moment of the Gaussian is compared on the
                    // scale of its standard deviation
                    const scalar difference =
                        mag(momentFromQuadrature - model.moments[mi])
                       /max(mag(model.moments[mi]), scalar(1e-3));

                    largest = max(largest, difference);
                }

                Info<< "  " << model.name << ", " << nMoments
                    << " moments, " << nNodes << " nodes: moments up to "
                    << "order " << 2*nNodes - 1
                    << " within " << largest << endl;

                if (largest > 1e-9)
                {
                    FatalErrorInFunction
                        << "GQMOM does not reproduce the " << model.name
                        << " distribution it continues the recurrence "
                        << "relationship with." << nl
                        << "    Number of moments: " << nMoments << nl
                        << "    Number of nodes: " << nNodes << nl
                        << "    Largest relative difference of a moment: "
                        << largest << nl
                        << exit(FatalError);
                }
            }
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

    inversion().invert(m);

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
//  The abscissae a quadrature fixes a node at scale with the coordinate as
//  every other abscissa does, so they are given here in the units of the
//  unscaled one, and written into the dictionary the quadrature is built
//  from at each scale. Gauss-Radau reads the lower one and Gauss-Lobatto
//  both; Gauss and GQMOM read neither.
void testScaleInvariance
(
    dictionary& dict,
    const word& quadratureName,
    const label nMaxNodes = 0,
    const scalar minKnownAbscissa = 0,
    const scalar maxKnownAbscissa = 1
)
{
    // Six moments of a quadrature of three nodes, which the inversion
    // reproduces exactly
    const scalarList w({0.2, 0.5, 0.3});
    const scalarList x({0.5, 1.5, 3.0});
    const label nMoments = 6;

    Info<< "\nTesting the scale invariance of " << quadratureName.c_str()
        << "\n" << endl;

    // The quadrature of the distribution written at one scale of its
    // abscissa, copied out of the inverter so that two can be compared
    auto quadratureAt = [&](const scalar scale)
    {
        dictionary scaledDict(dict);
        scaledDict.set("minKnownAbscissa", scale*minKnownAbscissa);
        scaledDict.set("maxKnownAbscissa", scale*maxKnownAbscissa);

        autoPtr<univariateMomentInversion> inversion
        (
            univariateMomentInversion::New(scaledDict, nMaxNodes)
        );

        const label nAdditionalPoints =
            inversion().nAdditionalQuadraturePoints(nMoments);

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

        inversion().invert(m);

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


// Gauss-Lobatto fixes a node at each of two abscissae, read from its
// dictionary. Given an interval that is empty, the system that corrects the
// last coefficients of the recurrence is singular, and the quadrature has to
// refuse the dictionary when it is built rather than return the 0/0 it used
// to, which surfaced later as a coefficient the Golub-Welsch algorithm found
// not to be finite.
//
// With exceptions on, any fatal error is caught, so the refusal is required
// to come from the Gauss-Lobatto quadrature itself.
void testLobattoRefusesAnEmptyInterval(dictionary& dict)
{
    Info<< "\nTesting that Gauss-Lobatto refuses an empty interval\n" << endl;

    // Both at zero; equal and away from zero; and the wrong way round
    const List<Pair<scalar>> intervals
    ({
        Pair<scalar>(0.0, 0.0),
        Pair<scalar>(1.5, 1.5),
        Pair<scalar>(3.5, 0.0)
    });

    for (const Pair<scalar>& interval : intervals)
    {
        dictionary lobattoDict(dict);
        lobattoDict.set("minKnownAbscissa", interval.first());
        lobattoDict.set("maxKnownAbscissa", interval.second());

        // The refusal is of a dictionary, which is an IO error
        const bool throwing = FatalError.throwing(true);
        const bool throwingIO = FatalIOError.throwing(true);

        bool refusedHere = false;
        string refusedElsewhere;

        try
        {
            univariateMomentInversion::New(lobattoDict);
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
        FatalIOError.throwing(throwingIO);

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
                << "building it failed in " << refusedElsewhere.c_str()
                << "." << nl
                << "    minKnownAbscissa: " << interval.first() << nl
                << "    maxKnownAbscissa: " << interval.second() << nl
                << exit(FatalError);
        }

        FatalErrorInFunction
            << "Gauss-Lobatto was built from a dictionary whose interval is "
            << "empty instead of refusing it." << nl
            << "    minKnownAbscissa: " << interval.first() << nl
            << "    maxKnownAbscissa: " << interval.second() << nl
            << exit(FatalError);
    }

    Info<< "\nGauss-Lobatto refuses an empty interval.\n" << endl;
}


// A moment set on the boundary of the moment space of [0, 1] has one measure,
// so every rule has to return it. A mixture fraction is such a set wherever
// a cell holds a pure stream, at zero or at one, beside anything else: two
// pure streams have m_k = m_1 for every k, and a Gauss rule used to give
// them one node at their mean, which is the two streams mixed, a Gauss-Radau
// rule invented fluid at zero where the pure stream was at one, and a
// Gauss-Lobatto rule stopped. The atoms are given, six moments are taken
// from them, and the inversion has to give the atoms back: those at both
// ends, at one end with interior ones, and interior ones only, which no
// boundary is involved in.
void testBoundaryOfTheUnitInterval(dictionary& dict, const string& name)
{
    Info<< "\nTesting " << name << " on the boundary of the moment space "
        << "of [0, 1]\n" << endl;

    const List<Tuple2<scalarList, scalarList>> measures
    ({
        {scalarList({0.7, 0.3}), scalarList({0, 1})},
        {scalarList({0.7, 0.3}), scalarList({0, 0.6})},
        {scalarList({0.7, 0.3}), scalarList({0.4, 1})},
        {scalarList({0.5, 0.3, 0.2}), scalarList({0, 0.5, 1})},
        {scalarList({0.5, 0.3, 0.2}), scalarList({0, 0.3, 0.8})},
        {scalarList({0.5, 0.3, 0.2}), scalarList({0.3, 0.8, 1})},
        {scalarList({0.7, 0.3}), scalarList({0.2, 0.7})}
    });

    for (const auto& measure : measures)
    {
        scalarList expectedWeights(measure.first());
        scalarList expectedAbscissae(measure.second());

        scalarList inputMoments(6, Zero);

        forAll(inputMoments, mi)
        {
            forAll(expectedWeights, nodei)
            {
                inputMoments[mi] +=
                    expectedWeights[nodei]*pow(expectedAbscissae[nodei], mi);
            }
        }

        univariateMomentSet m
        (
            inputMoments, supportType::ZeroOne, SMALL, SMALL
        );

        testQuadrature
        (
            m,
            dict,
            expectedWeights,
            expectedAbscissae,
            name + ", atoms at " + Foam::name(expectedAbscissae[0]) + " ... "
          + Foam::name(expectedAbscissae.last())
        );
    }
}


//- Two streams that are almost pure, a mixture fraction of 1e-6 and of
//  1 - 1e-6, have a second canonical moment within 4e-6 of one. With the
//  default boundaryTolerance they are not on the boundary of the moment
//  space at that order, and come back where they are. With a
//  boundaryTolerance of 1e-4 they are taken as pure: nodes at zero and at
//  one, with the weights that keep the moments of order zero and one.
void testBoundaryTolerance()
{
    Info<< "\nTesting the boundaryTolerance entry on [0, 1]\n" << endl;

    const scalar eps = 1e-6;
    const scalarList streamWeights({0.7, 0.3});
    const scalarList streamAbscissae({eps, 1 - eps});

    scalarList inputMoments(6, Zero);

    forAll(inputMoments, mi)
    {
        forAll(streamWeights, nodei)
        {
            inputMoments[mi] +=
                streamWeights[nodei]*pow(streamAbscissae[nodei], mi);
        }
    }

    for (const scalar tolerance : {scalar(-1), scalar(1e-4)})
    {
        const bool named = tolerance > 0;

        dictionary dict
        (
            IStringStream
            (
                "univariateMomentInversion Gauss;"
              + string
                (
                    named
                  ? " boundaryTolerance " + Foam::name(tolerance) + ";"
                  : ""
                )
            )()
        );

        autoPtr<univariateMomentInversion> inversion
        (
            univariateMomentInversion::New(dict)
        );

        univariateMomentSet m
        (
            inputMoments, supportType::ZeroOne, SMALL, SMALL
        );

        inversion().invert(m);

        const scalarList& weights(inversion().weights());
        const scalarList& abscissae(inversion().abscissae());

        Info<< "  boundaryTolerance " << inversion().boundaryTolerance()
            << ": weights " << weights << ", abscissae " << abscissae
            << endl;

        bool ok = (inversion().nNodes() == 2);

        if (ok && named)
        {
            // Pure streams, with the mass and the mean of the moments
            ok =
                mag(abscissae[0]) < 1e-14
             && mag(abscissae[1] - 1) < 1e-14
             && mag(weights[0] + weights[1] - inputMoments[0]) < 1e-14
             && mag(weights[1] - inputMoments[1]) < 1e-14;
        }
        else if (ok)
        {
            // The streams where they are. The inversion of two nodes this
            // close to the ends of the interval keeps about nine digits.
            ok =
                mag(abscissae[0] - eps) < 1e-9
             && mag(abscissae[1] - (1 - eps)) < 1e-9
             && mag(weights[0] - 0.7) < 1e-9
             && mag(weights[1] - 0.3) < 1e-9;
        }

        if (!ok)
        {
            FatalErrorInFunction
                << "The quadrature of two almost pure streams is not the "
                << "one the boundaryTolerance asks for." << nl
                << "    boundaryTolerance: "
                << inversion().boundaryTolerance() << nl
                << "    Weights: " << weights << nl
                << "    Abscissae: " << abscissae << nl
                << exit(FatalError);
        }
    }

    Info<< endl;
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


    // Test 6 - GQMOM on R+, ten moments of the uniform distribution on
    // [0, 1] continued with the gamma closure. There is no distribution to
    // compare with here: the values are those of the library, which an
    // independent implementation of the closure reproduces within 1e-13
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

    expectedWeightsTest6[0] = 0.03574921508537977;
    expectedWeightsTest6[1] = 0.09764650480516522;
    expectedWeightsTest6[2] = 0.1335134824107895;
    expectedWeightsTest6[3] = 0.1711684758894725;
    expectedWeightsTest6[4] = 0.1739860954488611;
    expectedWeightsTest6[5] = 0.1846183515715966;
    expectedWeightsTest6[6] = 0.1258528258568286;
    expectedWeightsTest6[7] = 0.07737017771086979;
    expectedWeightsTest6[8] = 9.485835768519413e-05;
    expectedWeightsTest6[9] = 1.286335127130966e-08;

    scalarList expectedAbscissaeTest6(10, 0);

    expectedAbscissaeTest6[0] = 0.01173191370526111;
    expectedAbscissaeTest6[1] = 0.08158521262637455;
    expectedAbscissaeTest6[2] = 0.19579690748922;
    expectedAbscissaeTest6[3] = 0.3519887349879595;
    expectedAbscissaeTest6[4] = 0.5232910050784757;
    expectedAbscissaeTest6[5] = 0.7063902877127978;
    expectedAbscissaeTest6[6] = 0.8652881472160519;
    expectedAbscissaeTest6[7] = 0.9662037783685198;
    expectedAbscissaeTest6[8] = 1.168841180426465;
    expectedAbscissaeTest6[9] = 1.557454260967021;

    testQuadrature
    (
        mGaussTest6,
        quadraturePropertiesGQMOM,
        expectedWeightsTest6,
        expectedAbscissaeTest6,
        "GQMOM",
        10
    );

    // Test 7 - GQMOM on [0, 1], ten moments of half a beta(3, 2)
    // distribution. The beta closure continues the canonical moments of a
    // beta distribution with its own, so the ten nodes have to be the
    // Gauss-Jacobi quadrature of the distribution, which is what the
    // expected values are: the roots of its orthogonal polynomial of degree
    // ten and their weights, found in 50-digit arithmetic. The quadrature
    // goes through the ten moments, whose round-off allows 1e-10.
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

    expectedWeightsTest7[0] = 0.00068702139102659208;
    expectedWeightsTest7[1] = 0.0064586770258096169;
    expectedWeightsTest7[2] = 0.024345252568145796;
    expectedWeightsTest7[3] = 0.056039175927690692;
    expectedWeightsTest7[4] = 0.09124606902902116;
    expectedWeightsTest7[5] = 0.11138498147045331;
    expectedWeightsTest7[6] = 0.10302860406764124;
    expectedWeightsTest7[7] = 0.069784170667645445;
    expectedWeightsTest7[8] = 0.030879685378606001;
    expectedWeightsTest7[9] = 0.0061463624739601487;

    scalarList expectedAbscissaeTest7(10, 0);

    expectedAbscissaeTest7[0] = 0.045250819584869484;
    expectedAbscissaeTest7[1] = 0.11844557799902641;
    expectedAbscissaeTest7[2] = 0.21738337124495932;
    expectedAbscissaeTest7[3] = 0.33538225149794593;
    expectedAbscissaeTest7[4] = 0.46441332586337908;
    expectedAbscissaeTest7[5] = 0.59568735822134655;
    expectedAbscissaeTest7[6] = 0.72025977492520121;
    expectedAbscissaeTest7[7] = 0.82964157493925428;
    expectedAbscissaeTest7[8] = 0.91637761443392429;
    expectedAbscissaeTest7[9] = 0.97454963563791952;

    testQuadrature
    (
        mGaussTest7,
        quadraturePropertiesGQMOM,
        expectedWeightsTest7,
        expectedAbscissaeTest7,
        "GQMOM",
        10,
        1e-10
    );

    // Test 8 - GQMOM on R+, the six moments case 8 of the aggregation and
    // breakup validation starts from, continued with the gamma closure to
    // ten nodes. As in test 6, the values are those of the library, which
    // an independent implementation of the closure reproduces within 1e-13
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

    expectedWeightsTest8[0] = 0.0001145883267614276;
    expectedWeightsTest8[1] = 0.0005787834611732504;
    expectedWeightsTest8[2] = 0.001257810942851072;
    expectedWeightsTest8[3] = 0.002468605049123289;
    expectedWeightsTest8[4] = 0.006133629917009222;
    expectedWeightsTest8[5] = 0.04525029932978637;
    expectedWeightsTest8[6] = 0.6772160600741152;
    expectedWeightsTest8[7] = 0.01487128917248305;
    expectedWeightsTest8[8] = 0.2521057336899993;
    expectedWeightsTest8[9] = 3.200036698790736e-06;

    scalarList expectedAbscissaeTest8(10, 0);

    expectedAbscissaeTest8[0] = 0.01622663939565731;
    expectedAbscissaeTest8[1] = 0.120727017195681;
    expectedAbscissaeTest8[2] = 0.2806738802943011;
    expectedAbscissaeTest8[3] = 0.4914640938694472;
    expectedAbscissaeTest8[4] = 0.7300805681269149;
    expectedAbscissaeTest8[5] = 0.971389438417085;
    expectedAbscissaeTest8[6] = 1.072903445254138;
    expectedAbscissaeTest8[7] = 1.231311738714221;
    expectedAbscissaeTest8[8] = 1.328905703996068;
    expectedAbscissaeTest8[9] = 1.534061102207788;

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
    testModelDistributionsGQMOM();

    testDegenerateVarianceGQMOM(quadraturePropertiesGQMOM, "gamma");
    testDegenerateVarianceGQMOM(quadraturePropertiesGQMOM, "lognormal");

    testScaleInvariance(quadraturePropertiesGauss, "Gauss");
    testScaleInvariance(quadraturePropertiesRadau, "Gauss-Radau");
    testScaleInvariance
    (
        quadraturePropertiesLobatto, "Gauss-Lobatto", 0, 0, 3.5
    );

    testLobattoRefusesAnEmptyInterval(quadraturePropertiesLobatto);
    testScaleInvariance(quadraturePropertiesGQMOM, "GQMOM", 5);

    testBoundaryTolerance();

    testBoundaryOfTheUnitInterval(quadraturePropertiesGauss, "Gauss");
    testBoundaryOfTheUnitInterval(quadraturePropertiesRadau, "Gauss-Radau");
    testBoundaryOfTheUnitInterval
    (
        quadraturePropertiesLobatto, "Gauss-Lobatto"
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
