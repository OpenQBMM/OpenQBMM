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
    Test-PopulationBalanceSubModels

Description
    Checks the source terms of the aggregation, breakup and growth models
    on a quadrature of three sizes, which reproduces its six moments
    exactly.

    Every expected value is computed here from the nodes, independently of
    the models: the fragments of each daughter distribution, the pairs that
    aggregate, and the rate of change of a moment under growth. The
    properties any correct model has are checked as well: breakup and
    aggregation conserve the volume of the particles, binary breakup adds
    one particle per event, aggregation with a constant kernel removes half
    of the pairs that meet, and growth changes a moment by its order times
    the moment below it.

    Run from the testCase directory, which holds a one-cell mesh, the
    moments, and the dictionaries of the models in subModelProperties.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "quadratureApproximations.H"
#include "aggregationKernel.H"
#include "breakupKernel.H"
#include "growthModel.H"
#include "daughterDistribution.H"

using namespace Foam;
using namespace Foam::populationBalanceSubModels;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nChecks = 0;
label nFailed = 0;


//- Compare a computed value with the one expected, relative to a scale
void compare
(
    const scalar computed,
    const scalar expected,
    const scalar scale,
    const string& what,
    const scalar tolerance = 1.0e-10
)
{
    nChecks++;

    const scalar error = mag(computed - expected)/max(mag(scale), VSMALL);

    // A value that is not a number fails the first comparison and an
    // infinite one the second
    if (computed == computed && mag(computed) <= VGREAT && error <= tolerance)
    {
        Info<< "  OK: \"" << what.c_str() << "\" (relative error "
            << error << ")" << nl;

        return;
    }

    nFailed++;

    Info<< "  FAILED: \"" << what.c_str() << "\"" << nl
        << "      computed " << computed << ", expected " << expected
        << ", relative error " << error << nl;
}


//- The moment of order k of the fragments a particle of size x breaks
//  into, derived from the fragments each distribution describes rather than
//  taken from its implementation
scalar fragments
(
    const word& distribution,
    const label k,
    const scalar x,
    const scalar primarySize
)
{
    // The global scope carries a second set of these, which makes a call
    // with plain scalars ambiguous
    using Foam::pow;
    using Foam::cbrt;
    using Foam::exp;

    if (distribution == "symmetricFragmentation")
    {
        // Two halves of the volume
        return 2.0*pow(x/cbrt(2.0), k);
    }
    else if (distribution == "uniform")
    {
        // Two fragments whose volume is uniform between zero and the parent
        return 6.0*pow(x, k)/(k + 3.0);
    }
    else if (distribution == "oneQuarterMassRatio")
    {
        // Four fifths and one fifth of the volume
        return (pow(0.8, k/3.0) + pow(0.2, k/3.0))*pow(x, k);
    }
    else if (distribution == "fullFragmentation")
    {
        // As many primary particles as the volume holds
        return pow3(x/primarySize)*pow(primarySize, k);
    }
    else if (distribution == "erosion")
    {
        // One primary particle and what remains, if there is anything to
        // erode
        if (x <= primarySize)
        {
            return pow(x, k);
        }

        return pow(primarySize, k) + pow(pow3(x) - pow3(primarySize), k/3.0);
    }

    FatalErrorInFunction
        << "No expected fragments for " << distribution << exit(FatalError);

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // The global scope carries a second set of these, which makes a call
    // with plain scalars ambiguous
    using Foam::pow;
    using Foam::cbrt;
    using Foam::exp;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    scalarQuadratureApproximation quadrature
    (
        "populationBalance",
        mesh,
        List<supportType>(1, supportType::RPlus)
    );

    quadrature.updateQuadrature();

    IOdictionary subModels
    (
        IOobject
        (
            "subModelProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const label celli = 0;

    // The nodes every model sees, in increasing size
    const label nNodes = quadrature.nodes().size();
    scalarList w(nNodes);
    scalarList x(nNodes);

    forAll(w, nodei)
    {
        w[nodei] = quadrature.nodes()[nodei].weight()[celli];
        x[nodei] = quadrature.nodes()[nodei].abscissae()[0][celli];
    }

    auto moment = [&](const scalar k)
    {
        scalar m = 0;

        forAll(w, i)
        {
            m += w[i]*pow(x[i], k);
        }

        return m;
    };

    Info<< "\nThe quadrature" << endl;
    {
        labelList order;
        sortedOrder(x, order);

        const scalarList xCase({1.0e-5, 2.0e-5, 4.0e-5});
        const scalarList wCase({2.0e11, 5.0e11, 3.0e11});

        forAll(order, i)
        {
            compare
            (
                x[order[i]], xCase[i], xCase[i],
                "size " + Foam::name(i) + " is the one the case was built of",
                1.0e-8
            );

            compare
            (
                w[order[i]], wCase[i], wCase[i],
                "weight " + Foam::name(i) + " is the one the case was built of",
                1.0e-8
            );
        }
    }


    // * * * * * * * * * * * * * * * Aggregation * * * * * * * * * * * * * //

    Info<< "\nAggregation" << endl;

    for (const word kernelName : {"constantAggregation", "sumAggregation"})
    {
        const dictionary& dict = subModels.subDict(kernelName);
        autoPtr<aggregationKernel> kernel(aggregationKernel::New(dict, mesh));

        const scalar Ca = dimensionedScalar("Ca", dict).value();
        const bool sum = (kernelName == "sumAggregation");

        auto Ka = [&](const scalar xi, const scalar xj)
        {
            return sum ? Ca*(pow3(xi) + pow3(xj)) : Ca;
        };

        for (label k = 0; k <= 5; k++)
        {
            scalar birth = 0;
            scalar death = 0;

            forAll(w, i)
            {
                forAll(w, j)
                {
                    birth +=
                        w[i]*w[j]*Ka(x[i], x[j])
                       *0.5*pow(pow3(x[i]) + pow3(x[j]), k/3.0);

                    death += w[i]*w[j]*Ka(x[i], x[j])*pow(x[i], k);
                }
            }

            compare
            (
                kernel->aggregationSource
                (
                    labelList(1, k), celli, quadrature, 0
                ),
                birth - death,
                max(birth, death),
                kernelName + ", moment of order " + Foam::name(k)
            );
        }

        const scalar m0 = moment(0);
        const scalar m3 = moment(3);

        compare
        (
            kernel->aggregationSource(labelList(1, 3), celli, quadrature, 0),
            0,
            (sum ? Ca*sqr(m3) : Ca*m0*m3),
            kernelName + " conserves the volume of the particles"
        );

        compare
        (
            kernel->aggregationSource(labelList(1, 0), celli, quadrature, 0),
            (sum ? -Ca*m0*m3 : -0.5*Ca*sqr(m0)),
            (sum ? Ca*m0*m3 : 0.5*Ca*sqr(m0)),
            kernelName + " removes the particles of the pairs that meet"
        );
    }


    // * * * * * * * * * * * * * * * * Breakup  * * * * * * * * * * * * * * //

    Info<< "\nBreakup" << endl;

    // Each daughter distribution, with a constant kernel and no smallest
    // size, so that every size breaks
    for (const entry& e : subModels.subDict("breakup"))
    {
        const dictionary& dict = e.dict();
        autoPtr<breakupKernel> kernel(breakupKernel::New(dict, mesh));

        const dictionary& daughterDict = dict.subDict("daughterDistribution");
        const word distribution(daughterDict.get<word>("daughterDistribution"));

        const scalar primarySize =
            daughterDict.found("primarySize")
          ? dimensionedScalar("primarySize", daughterDict).value()
          : scalar(0);

        const scalar Cb = dimensionedScalar("Cb", dict).value();

        for (label k = 0; k <= 5; k++)
        {
            scalar expected = 0;
            scalar scale = 0;

            forAll(w, i)
            {
                expected +=
                    w[i]*Cb
                   *(
                        fragments(distribution, k, x[i], primarySize)
                      - pow(x[i], k)
                    );

                scale += w[i]*Cb*pow(x[i], k);
            }

            compare
            (
                kernel->breakupSource(labelList(1, k), celli, quadrature),
                expected,
                scale,
                distribution + ", moment of order " + Foam::name(k)
            );
        }

        compare
        (
            kernel->breakupSource(labelList(1, 3), celli, quadrature),
            0,
            Cb*moment(3),
            distribution + " conserves the volume of the particles"
        );

        if (fragments(distribution, 0, x[0], primarySize) == 2)
        {
            compare
            (
                kernel->breakupSource(labelList(1, 0), celli, quadrature),
                Cb*moment(0),
                Cb*moment(0),
                distribution + " adds one particle per event"
            );
        }

        // The distribution alone, on either coordinate, including a size
        // below that of a primary particle
        autoPtr<daughterDistribution> daughters
        (
            daughterDistribution::New(daughterDict)
        );

        for (const scalar size : {2.0e-6, 1.0e-5, 4.0e-5})
        {
            compare
            (
                daughters->mD(3, size),
                pow3(size),
                pow3(size),
                distribution + " conserves volume on a length coordinate, "
              + "size " + Foam::name(size)
            );

            compare
            (
                daughters->mDMass(1, size),
                size,
                size,
                distribution + " conserves mass on a mass coordinate, "
              + "size " + Foam::name(size)
            );
        }
    }

    // A smallest size between the nodes, then the other two kernels
    {
        const dictionary& dict = subModels.subDict("gatedBreakup");
        autoPtr<breakupKernel> kernel(breakupKernel::New(dict, mesh));

        const scalar Cb = dimensionedScalar("Cb", dict).value();
        const scalar minAbscissa = dict.get<scalar>("minAbscissa");

        for (label k = 0; k <= 5; k++)
        {
            scalar expected = 0;
            scalar scale = 0;

            forAll(w, i)
            {
                if (x[i] >= minAbscissa)
                {
                    expected +=
                        w[i]*Cb
                       *(fragments("symmetricFragmentation", k, x[i], 0)
                       - pow(x[i], k));
                }

                scale += w[i]*Cb*pow(x[i], k);
            }

            compare
            (
                kernel->breakupSource(labelList(1, k), celli, quadrature),
                expected,
                scale,
                "only the sizes above minAbscissa break, moment of order "
              + Foam::name(k)
            );
        }
    }

    for (const word kernelName : {"powerLawBreakup", "exponentialBreakup"})
    {
        const dictionary& dict = subModels.subDict(kernelName);
        autoPtr<breakupKernel> kernel(breakupKernel::New(dict, mesh));

        const scalar Cb = dimensionedScalar("Cb", dict).value();
        const bool powerLaw = (kernelName == "powerLawBreakup");

        auto Kb = [&](const scalar size)
        {
            if (powerLaw)
            {
                return Cb*pow(size, dict.get<scalar>("abscissaExponent"));
            }

            return
                Cb*exp(dimensionedScalar("expCoeff", dict).value()*pow3(size));
        };

        for (label k = 0; k <= 5; k++)
        {
            scalar expected = 0;
            scalar scale = 0;

            forAll(w, i)
            {
                expected +=
                    w[i]*Kb(x[i])
                   *(fragments("symmetricFragmentation", k, x[i], 0)
                   - pow(x[i], k));

                scale += w[i]*Kb(x[i])*pow(x[i], k);
            }

            compare
            (
                kernel->breakupSource(labelList(1, k), celli, quadrature),
                expected,
                scale,
                kernelName + ", moment of order " + Foam::name(k)
            );
        }
    }


    // * * * * * * * * * * * * * * * * Growth * * * * * * * * * * * * * * * //

    Info<< "\nGrowth" << endl;

    {
        const dictionary& dict = subModels.subDict("constantGrowth");
        autoPtr<growthModel> model(growthModel::New(dict, mesh));

        const scalar Cg = dict.get<scalar>("Cg");

        compare
        (
            model->phaseSpaceConvection(labelList(1, 0), celli, quadrature),
            0,
            Cg*moment(0),
            "constant growth leaves the number of particles unchanged"
        );

        for (label k = 1; k <= 5; k++)
        {
            compare
            (
                model->phaseSpaceConvection(labelList(1, k), celli, quadrature),
                k*Cg*moment(k - 1),
                k*Cg*moment(k - 1),
                "constant growth changes the moment of order "
              + Foam::name(k) + " by its order times the one below"
            );
        }
    }

    // The bounds, on the model that always had them and on one that did not
    for (const word modelName : {"gatedGrowth", "boundedEvaporation"})
    {
        const dictionary& dict = subModels.subDict(modelName);
        autoPtr<growthModel> model(growthModel::New(dict, mesh));

        const scalar Cg = dict.get<scalar>("Cg");
        const scalar minAbscissa = dict.getOrDefault<scalar>("minAbscissa", 0);
        const scalar maxAbscissa =
            dict.getOrDefault<scalar>("maxAbscissa", GREAT);

        const bool evaporation = (modelName == "boundedEvaporation");

        // The linear evaporation rate of a length, from the rate of change
        // of the volume of a sphere
        auto Kg = [&](const scalar d)
        {
            if (evaporation)
            {
                return
                   -Cg*constant::mathematical::pi/6.0*pow3(d)
                   *2.0/(constant::mathematical::pi*sqr(d));
            }

            return Cg;
        };

        for (label k = 1; k <= 5; k++)
        {
            scalar expected = 0;
            scalar scale = 0;

            forAll(w, i)
            {
                const scalar term = w[i]*Kg(x[i])*k*pow(x[i], k - 1);

                if (x[i] >= minAbscissa && x[i] <= maxAbscissa)
                {
                    expected += term;
                }

                scale += mag(term);
            }

            compare
            (
                model->phaseSpaceConvection(labelList(1, k), celli, quadrature),
                expected,
                scale,
                modelName + " grows only the sizes within its bounds, "
              + "moment of order " + Foam::name(k)
            );
        }
    }

    // Bounds and a coefficient written the older way, with a name and
    // dimensions, as thirty of the cases write them
    {
        const dictionary& dict = subModels.subDict("legacyGrowth");
        autoPtr<growthModel> model(growthModel::New(dict, mesh));

        compare
        (
            model->phaseSpaceConvection(labelList(1, 1), celli, quadrature),
            moment(0),
            moment(0),
            "a dictionary written with names and dimensions is read"
        );
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;

    if (nFailed > 0)
    {
        FatalErrorInFunction
            << nFailed << " of " << nChecks << " checks failed."
            << exit(FatalError);
    }

    Info<< nChecks << " checks passed." << nl << endl;

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
