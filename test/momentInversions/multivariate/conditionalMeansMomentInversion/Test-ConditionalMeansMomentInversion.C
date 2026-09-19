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
    Test-ConditionalMeansMomentInversion

Description
    Test the inversion of a distribution in one primary coordinate and
    coordinates conditioned on it, each carried by its mean at every node.

    A distribution of that form, n nodes of the primary coordinate each
    carrying one value of every other coordinate, is its own quadrature, so
    the inversion has to recover its nodes, and every moment of it, to
    round-off. Each check builds one from known nodes, with the moment set
    the inversion builds from counts, and inverts it:

    - a mixture fraction with two compositions, every coordinate on [0, 1];
    - a primary coordinate on the reals and on the positive reals;
    - a size with a conditioned scalar and a conditioned velocity;
    - a node carrying less than smallM0, which is not read back while the
      others are still recovered exactly;
    - conditioned values outside [0, 1] or below zero, which are clipped to
      the support of their dimension, and counted, and left alone on the
      reals;
    - the selection of the inversion by both of its names.

    The test is self-contained and reads nothing from the disk.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IStringStream.H"
#include "conditionalMeansMomentInversion.H"
#include "multivariateMomentSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nChecks = 0;
label nFailed = 0;

void check(const word& name, const bool ok, const string& detail = "")
{
    nChecks++;

    if (ok)
    {
        Info<< "    passed: " << name << nl;
    }
    else
    {
        nFailed++;
        Info<< "    FAILED: " << name << " " << detail.c_str() << nl;
    }
}


bool near(const scalar a, const scalar b, const scalar tol = 1e-10)
{
    return mag(a - b) <= tol*max(1.0, max(mag(a), mag(b)));
}


dictionary inversionDict(const word& type = "conditionalMeans")
{
    return dictionary
    (
        IStringStream
        (
            "type " + type + ";"
            "basicQuadrature { univariateMomentInversion Gauss; }"
        )()
    );
}


//- A distribution of nodes of the primary coordinate xs with weights ws,
//  each carrying the value phis[d][node] of conditioned dimension d + 1
struct Distribution
{
    scalarList ws;
    scalarList xs;
    List<scalarList> phis;
};


//- The moments of a distribution for the given orders
multivariateMomentSet momentsOf
(
    const Distribution& dist,
    const labelListList& orders,
    const List<supportType>& supports
)
{
    multivariateMomentSet moments(orders.size(), orders, supports, SMALL, 0);

    forAll(orders, mi)
    {
        const labelList& order = orders[mi];
        scalar m = 0;

        forAll(dist.ws, nodei)
        {
            scalar c = dist.ws[nodei]*pow(dist.xs[nodei], order[0]);

            forAll(dist.phis, d)
            {
                c *= pow(dist.phis[d][nodei], order[d + 1]);
            }

            m += c;
        }

        moments(order) = m;
    }

    return moments;
}


//- Invert a distribution and compare what comes out with what is
//  expected: the weights, the primary abscissae and the conditioned values,
//  the velocity ones among them given by velocityIndexes
void recover
(
    const word& name,
    const Distribution& dist,
    const List<supportType>& supports,
    const label nScalars,
    const labelList& velocityIndexes,
    const List<scalarList>& expected,
    const label expectedClipped = 0
)
{
    Info<< nl << name << nl;

    const label nVelocities =
        (velocityIndexes.size() == 1 && velocityIndexes[0] < 0)
      ? 0
      : velocityIndexes.size();

    labelListList orders;
    labelListList nodes;

    multivariateMomentInversions::conditionalMeans::momentOrdersAndNodes
    (
        2*dist.ws.size(),
        nScalars,
        nVelocities,
        orders,
        nodes
    );

    const multivariateMomentSet moments(momentsOf(dist, orders, supports));

    multivariateMomentInversions::conditionalMeans inverter
    (
        inversionDict(), orders, nodes, velocityIndexes
    );

    check("the inversion succeeds", inverter.invert(moments));

    const mappedScalarList& weights = inverter.weights();
    const mappedList<scalarList>& abscissae = inverter.abscissae();
    const mappedVectorList& velocities = inverter.velocityAbscissae();

    bool weightsOk = true;
    bool primaryOk = true;
    bool conditionedOk = true;

    forAll(nodes, nodei)
    {
        // The weights are compared relative to the largest, since those of
        // a distribution can be of any scale
        weightsOk =
            weightsOk
         && mag(weights[nodei] - dist.ws[nodei]) <= 1e-10*max(dist.ws);

        primaryOk = primaryOk && near(abscissae[nodei][0], dist.xs[nodei]);

        label scalari = 1;
        label velocityi = 0;

        forAll(expected, d)
        {
            const label dimi = d + 1;
            const bool isVelocity = velocityIndexes.found(dimi);

            const scalar value =
                isVelocity
              ? velocities[nodei][velocityi++]
              : abscissae[nodei][scalari++];

            if (!near(value, expected[d][nodei]))
            {
                conditionedOk = false;

                Info<< "        node " << nodei << ", dimension " << dimi
                    << ": " << value << ", expected " << expected[d][nodei]
                    << nl;
            }
        }
    }

    check("the weights", weightsOk);
    check("the primary abscissae", primaryOk);
    check("the conditioned means", conditionedOk);
    check
    (
        "the number of clipped values",
        inverter.nClipped() == expectedClipped,
        "(" + Foam::name(inverter.nClipped()) + ", expected "
      + Foam::name(expectedClipped) + ")"
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"

    const supportType R = supportType::R;
    const supportType RPlus = supportType::RPlus;
    const supportType ZeroOne = supportType::ZeroOne;
    const labelList noVelocity({-1});

    // A mixture fraction with two compositions
    {
        Distribution dist;
        dist.ws = scalarList({0.2, 0.5, 0.3});
        dist.xs = scalarList({0.1, 0.45, 0.9});
        dist.phis = List<scalarList>
        ({
            scalarList({0.8, 0.3, 0.05}),
            scalarList({0.0, 0.6, 1.0})
        });

        recover
        (
            "A mixture fraction with two compositions, all on [0, 1]",
            dist, List<supportType>({ZeroOne, ZeroOne, ZeroOne}), 2,
            noVelocity, dist.phis
        );
    }

    // The primary coordinate on the reals and on the positive reals
    {
        Distribution dist;
        dist.ws = scalarList({0.3, 0.3, 0.4});
        dist.xs = scalarList({-1.2, 0.3, 2.0});
        dist.phis = List<scalarList>({scalarList({0.25, 0.5, 0.75})});

        recover
        (
            "A primary coordinate on the reals",
            dist, List<supportType>({R, ZeroOne}), 1, noVelocity, dist.phis
        );

        dist.xs = scalarList({0.5, 2.0, 7.0});

        recover
        (
            "A primary coordinate on the positive reals",
            dist, List<supportType>({RPlus, ZeroOne}), 1, noVelocity,
            dist.phis
        );
    }

    // A size, with a scalar and a velocity conditioned on it: the scalar
    // is stored with the abscissae and the velocity with the velocities
    {
        Distribution dist;
        dist.ws = scalarList({1e6, 2e6, 5e5});
        dist.xs = scalarList({1e-4, 3e-4, 8e-4});
        dist.phis = List<scalarList>
        ({
            scalarList({0.9, 0.4, 0.1}),
            scalarList({-0.3, 0.2, 1.5})
        });

        recover
        (
            "A size with a conditioned scalar and velocity",
            dist, List<supportType>({RPlus, ZeroOne, R}), 1,
            labelList({2}), dist.phis
        );
    }

    // A node that carries less than smallM0 (SMALL, absolute) stays in the
    // system, so that the others are recovered exactly, and is not read
    // back. Its weight is small against smallM0 but not against the other
    // nodes, so that the moments still resolve it.
    {
        Distribution dist;
        dist.ws = scalarList({1e-14, 5e-16, 1e-14});
        dist.xs = scalarList({0.2, 0.5, 0.8});
        dist.phis = List<scalarList>({scalarList({0.3, 0.9, 0.6})});

        recover
        (
            "A node carrying less than smallM0",
            dist, List<supportType>({ZeroOne, ZeroOne}), 1, noVelocity,
            List<scalarList>({scalarList({0.3, 0.0, 0.6})})
        );
    }

    // Conditioned values outside their support are clipped to it, and
    // counted, and left alone on the reals
    {
        Distribution dist;
        dist.ws = scalarList({0.3, 0.4, 0.3});
        dist.xs = scalarList({0.2, 0.5, 0.8});
        dist.phis = List<scalarList>({scalarList({1.2, 0.5, -0.1})});

        recover
        (
            "Values outside [0, 1] on [0, 1]",
            dist, List<supportType>({ZeroOne, ZeroOne}), 1, noVelocity,
            List<scalarList>({scalarList({1.0, 0.5, 0.0})}), 2
        );

        recover
        (
            "The same values on the positive reals",
            dist, List<supportType>({ZeroOne, RPlus}), 1, noVelocity,
            List<scalarList>({scalarList({1.2, 0.5, 0.0})}), 1
        );

        recover
        (
            "The same values on the reals",
            dist, List<supportType>({ZeroOne, R}), 1, noVelocity,
            dist.phis, 0
        );
    }

    // Both names select the same inversion
    {
        Info<< nl << "Selection by name" << nl;

        labelListList orders;
        labelListList nodes;

        multivariateMomentInversions::conditionalMeans::momentOrdersAndNodes
        (
            4, 0, 1, orders, nodes
        );

        forAll(wordList({"conditionalMeans", "monoKinetic"}), i)
        {
            const word type(wordList({"conditionalMeans", "monoKinetic"})[i]);

            autoPtr<multivariateMomentInversion> inverter
            (
                multivariateMomentInversion::New
                (
                    inversionDict(type), orders, nodes, labelList({1})
                )
            );

            check
            (
                type + " selects conditionalMeans",
                inverter->type() == "conditionalMeans"
            );
        }
    }

    Info<< nl;

    if (nFailed > 0)
    {
        FatalErrorInFunction
            << nFailed << " of " << nChecks << " checks failed." << nl
            << exit(FatalError);
    }

    Info<< "All " << nChecks << " checks passed" << nl << endl;

    return 0;
}

// ************************************************************************* //
