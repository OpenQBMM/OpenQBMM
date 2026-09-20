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
    along with OpenQBMM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-MixtureFractionStreams

Description
    Test the moments generated for a mixture fraction and the compositions
    conditioned on it, where two streams have only mixed.

    Three species, A and B in the stream at xi = 0 and C alone in the one
    at xi = 1, and the moment set of two nodes: <xi^k> for k < 4 and
    <xi^k Y_j> for k < 2. For each distribution of the mixture fraction:
    - beta: its moments are those of the beta function, written with the
      gamma function and not with the recurrence the model uses, its mean
      and its variance are the ones asked for, <xi^k Y_C> is <xi^(k+1)>
      since C comes only with the stream at one, <xi^k Y_A> is a quarter
      of <xi^k> - <xi^(k+1)>, and the compositions sum to <xi^k>;
    - delta: at zero the cell holds the first stream and at one the second;
    - nodes: every moment is the sum over the nodes of the weight, the
      power of the abscissa and the composition at that abscissa, which
      does not go through the moments of the mixture fraction.
    With one coordinate more the distribution carries the enthalpy, and
    its moments are those of an enthalpy linear in the mixture fraction
    between the streams, whose enthalpies the model takes from the thermo
    of the case: they are checked against the heats of formation and the
    constant heat capacities of the species, written out here.
    The moments of a patch and of a set of cells have the size of what
    they are for. A moment of second order in a composition, and a variance
    a beta distribution cannot have, are refused.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hexCellFvMesh.H"
#include "momentGenerationModel.H"
#include "StringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nChecks = 0;
label nFailed = 0;

void check(const word& name, const scalar error, const scalar tolerance)
{
    nChecks++;

    if (error <= tolerance)
    {
        Info<< "    passed: " << name << " (" << error << ")" << nl;
    }
    else
    {
        nFailed++;
        Info<< "    FAILED: " << name << " (" << error << ", tolerance "
            << tolerance << ")" << nl;
    }
}


labelList order(const label k, const label composition = -1)
{
    labelList o(4, 0);
    o[0] = k;

    if (composition >= 0)
    {
        o[1 + composition] = 1;
    }

    return o;
}


int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"

    Foam::simplifiedMeshes::hexCellFvMesh mesh(runTime);

    // The moment set of two nodes of the mixture fraction and three
    // compositions
    DynamicList<labelList> orders;

    for (label k = 0; k < 4; k++)
    {
        orders.append(order(k));
    }

    for (label j = 0; j < 3; j++)
    {
        for (label k = 0; k < 2; k++)
        {
            orders.append(order(k, j));
        }
    }

    const labelListList momentOrders(orders);

    const string modelText
    (
        "type mixtureFractionStreams;"
        "weightDimension [0 0 0 0 0 0 0];"
        "abscissaDimension [0 0 0 0 0 0 0];"
        "species (A B C);"
        "streams"
        "{"
        "    xi0 { T 300; Y { A 0.25; B 0.75; } }"
        "    xi1 { T 300; Y { C 1; } }"
        "}"
    );

    IStringStream modelStream(modelText);
    const dictionary modelDict(modelStream);

    autoPtr<momentGenerationModel> generator
    (
        momentGenerationModel::New(mesh, modelDict, "mixing", momentOrders, 2)
    );

    const scalarField Y0({0.25, 0.75, 0});
    const scalarField Y1({0, 0, 1});

    auto moment = [&generator](const labelList& o)
    {
        return generator->moments()(o)[0];
    };


    // * * * * * * * * * * * * * * * * * Beta  * * * * * * * * * * * * * * //

    Info<< nl << "Beta distribution of mean 0.3 and variance 0.05" << nl;
    {
        const scalar mean = 0.3;
        const scalar variance = 0.05;

        generator->updateMoments
        (
            dictionary
            (
                IStringStream
                (
                    "distribution beta; mean 0.3; variance 0.05;"
                )()
            )
        );

        const scalar c = mean*(1 - mean)/variance - 1;
        const scalar a = mean*c;
        const scalar b = (1 - mean)*c;

        scalar worstBeta = 0;

        for (label k = 0; k < 4; k++)
        {
            const scalar exact =
                Foam::exp
                (
                    Foam::lgamma(a + k) + Foam::lgamma(a + b)
                  - Foam::lgamma(a) - Foam::lgamma(a + b + k)
                );

            worstBeta = max(worstBeta, mag(moment(order(k)) - exact));
        }

        check("the moments are those of the beta function", worstBeta, 1e-13);

        check
        (
            "the mean is the one asked for",
            mag(moment(order(1)) - mean),
            1e-14
        );

        check
        (
            "the variance is the one asked for",
            mag(moment(order(2)) - sqr(moment(order(1))) - variance),
            1e-14
        );

        scalar worstC = 0;
        scalar worstA = 0;
        scalar worstSum = 0;

        for (label k = 0; k < 2; k++)
        {
            const scalar mk = moment(order(k));
            const scalar mk1 = moment(order(k + 1));

            worstC = max(worstC, mag(moment(order(k, 2)) - mk1));
            worstA = max(worstA, mag(moment(order(k, 0)) - 0.25*(mk - mk1)));

            worstSum =
                max
                (
                    worstSum,
                    mag
                    (
                        moment(order(k, 0)) + moment(order(k, 1))
                      + moment(order(k, 2)) - mk
                    )
                );
        }

        check("<xi^k Y_C> is <xi^(k+1)>", worstC, 1e-15);
        check("<xi^k Y_A> is a quarter of <xi^k> - <xi^(k+1)>", worstA, 1e-15);
        check("the compositions sum to <xi^k>", worstSum, 1e-15);
    }


    // * * * * * * * * * * * * * * * * * Delta * * * * * * * * * * * * * * //

    Info<< nl << "All of the fluid at one mixture fraction" << nl;
    {
        generator->updateMoments
        (
            dictionary(IStringStream("distribution delta; value 0;")())
        );

        scalar worst = mag(moment(order(0)) - 1);

        for (label k = 1; k < 4; k++)
        {
            worst = max(worst, mag(moment(order(k))));
        }

        for (label j = 0; j < 3; j++)
        {
            worst = max(worst, mag(moment(order(0, j)) - Y0[j]));
            worst = max(worst, mag(moment(order(1, j))));
        }

        check("at zero the cell holds the first stream", worst, 0);

        generator->updateMoments
        (
            dictionary(IStringStream("distribution delta; value 1;")())
        );

        worst = 0;

        for (label k = 0; k < 4; k++)
        {
            worst = max(worst, mag(moment(order(k)) - 1));
        }

        for (label j = 0; j < 3; j++)
        {
            worst = max(worst, mag(moment(order(0, j)) - Y1[j]));
            worst = max(worst, mag(moment(order(1, j)) - Y1[j]));
        }

        check("at one the cell holds the second stream", worst, 0);
    }


    // * * * * * * * * * * * * * * * * * Nodes * * * * * * * * * * * * * * //

    Info<< nl << "Two nodes of the mixture fraction" << nl;
    {
        const scalarList w({0.4, 0.6});
        const scalarList xi({0.2, 0.9});

        IStringStream nodesStream
        (
            "distribution nodes; weights (0.4 0.6); abscissae (0.2 0.9);"
        );
        const dictionary nodesDict(nodesStream);

        generator->updateMoments(nodesDict);

        scalar worst = 0;

        forAll(momentOrders, mi)
        {
            const labelList& o = momentOrders[mi];
            scalar exact = 0;

            forAll(w, nodei)
            {
                scalar term = w[nodei]*pow(xi[nodei], o[0]);

                for (label j = 0; j < 3; j++)
                {
                    if (o[1 + j] == 1)
                    {
                        term *= xi[nodei]*Y1[j] + (1 - xi[nodei])*Y0[j];
                    }
                }

                exact += term;
            }

            worst = max(worst, mag(moment(o) - exact));
        }

        check("every moment is the sum over the nodes", worst, 1e-15);

        // The sizes of a patch and of a set of cells
        generator->updateMoments(nodesDict, 0);

        check
        (
            "the moments of a patch have the size of the patch",
            mag
            (
                generator->moments()(order(0)).size()
              - mesh.boundaryMesh()[0].size()
            ),
            0
        );

        generator->updateMoments(nodesDict, labelList({0}));

        check
        (
            "the moments of a set of cells have the size of the set",
            mag(generator->moments()(order(1)).size() - 1)
          + mag(moment(order(1)) - (0.4*0.2 + 0.6*0.9)),
            1e-15
        );
    }


    // * * * * * * * * * * * * * * * * Enthalpy  * * * * * * * * * * * * * //

    Info<< nl << "The enthalpy conditioned on the mixture fraction" << nl;
    {
        // The same set with one coordinate more, the enthalpy
        DynamicList<labelList> hOrders;

        for (label k = 0; k < 4; k++)
        {
            labelList o(5, 0);
            o[0] = k;
            hOrders.append(o);
        }

        for (label j = 0; j < 4; j++)
        {
            for (label k = 0; k < 2; k++)
            {
                labelList o(5, 0);
                o[0] = k;
                o[1 + j] = 1;
                hOrders.append(o);
            }
        }

        const string hModelText
        (
            "type mixtureFractionStreams;"
            "weightDimension [0 0 0 0 0 0 0];"
            "abscissaDimension [0 0 0 0 0 0 0];"
            "species (A B C);"
            "streams"
            "{"
            "    xi0 { T 400; Y { A 0.25; B 0.75; } }"
            "    xi1 { T 300; Y { C 1; } }"
            "}"
        );

        IStringStream hModelStream(hModelText);
        const dictionary hModelDict(hModelStream);

        autoPtr<momentGenerationModel> hGenerator
        (
            momentGenerationModel::New
            (
                mesh,
                hModelDict,
                "mixing",
                labelListList(hOrders),
                2
            )
        );

        hGenerator->updateMoments
        (
            dictionary
            (
                IStringStream
                (
                    "distribution beta; mean 0.3; variance 0.05;"
                )()
            )
        );

        // The absolute enthalpy of a species of constant heat capacity is
        // its heat of formation and Cp (T - Tstd)
        const scalar Tstd = 298.15;
        const scalar hA = 0 + 1000*(400 - Tstd);
        const scalar hB = 100000 + 1100*(400 - Tstd);
        const scalar hC = -200000 + 1200*(300 - Tstd);
        const scalar h0 = 0.25*hA + 0.75*hB;
        const scalar h1 = hC;

        scalar worst = 0;

        for (label k = 0; k < 2; k++)
        {
            labelList xiOrder(5, 0);
            labelList hOrder(5, 0);
            xiOrder[0] = k;
            hOrder[0] = k;
            hOrder[4] = 1;

            const scalar mk = hGenerator->moments()(xiOrder)[0];
            xiOrder[0] = k + 1;
            const scalar mk1 = hGenerator->moments()(xiOrder)[0];

            worst =
                max
                (
                    worst,
                    mag
                    (
                        hGenerator->moments()(hOrder)[0]
                      - (h1*mk1 + h0*(mk - mk1))
                    )/mag(h1 - h0)
                );
        }

        check("<xi^k h> is linear between the streams", worst, 1e-14);

        labelList hOrder(5, 0);
        hOrder[4] = 1;

        check
        (
            "a moment of the enthalpy is an energy per unit mass",
            hGenerator->momentDims()(hOrder) == dimEnergy/dimMass ? 0 : 1,
            0
        );
    }


    // * * * * * * * * * * * * * * * * Refusals  * * * * * * * * * * * * * //

    Info<< nl << "What the model refuses" << nl;
    {
        const bool oldThrowing = FatalIOError.throwing(true);

        bool refusedVariance = false;

        try
        {
            generator->updateMoments
            (
                dictionary
                (
                    IStringStream
                    (
                        "distribution beta; mean 0.3; variance 0.3;"
                    )()
                )
            );
        }
        catch (const Foam::IOerror&)
        {
            refusedVariance = true;
        }

        check
        (
            "a variance a beta distribution cannot have",
            refusedVariance ? 0 : 1,
            0
        );

        bool refusedOrder = false;

        try
        {
            labelListList badOrders(momentOrders);
            badOrders.last()[3] = 2;

            autoPtr<momentGenerationModel> bad
            (
                momentGenerationModel::New
                (
                    mesh,
                    modelDict,
                    "mixing",
                    badOrders,
                    2
                )
            );
        }
        catch (const Foam::IOerror&)
        {
            refusedOrder = true;
        }

        check
        (
            "a moment of second order in a composition",
            refusedOrder ? 0 : 1,
            0
        );

        FatalIOError.throwing(oldThrowing);
    }

    Info<< nl;

    if (nFailed)
    {
        Info<< nFailed << " of " << nChecks << " checks FAILED" << nl << endl;
        return 1;
    }

    Info<< "All " << nChecks << " checks passed" << nl << endl;

    return 0;
}


// ************************************************************************* //
