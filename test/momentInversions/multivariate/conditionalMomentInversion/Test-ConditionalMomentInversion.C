/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2026 Alberto Passalacqua
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
    Test-ConditionalMomentInversion

Description
    Test the conditional quadrature method of moments.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "mappedLists.H"
#include "supportType.H"
#include "conditionalMomentInversion.H"
#include "Random.H"
#include <cmath>
#include "multivariateMomentTest.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A conditional quadrature, which is what the method represents exactly.
//
//  The weights and abscissae of a direction depend on the node the
//  preceding directions are at, and on nothing else. A distribution built
//  that way is its own conditional quadrature, so the inversion recovers it
//  and every moment of it has to come back.
//
//  A quadrature drawn at random is not of that form: its abscissae in a
//  direction depend on the whole node index rather than on the preceding
//  part of it, so the method closes it rather than reproducing it, and
//  moments come back with tens of percent of error while it behaves as
//  designed.

//- A number that depends on the node the directions before dimi are at
scalar prefixKey(const labelList& nodeIndex, const label dimi)
{
    scalar key = 0.0;

    for (label i = 0; i < dimi; i++)
    {
        key = 3.0*key + nodeIndex[i];
    }

    return key;
}


//- Weight of the node of direction dimi, conditioned on the ones before it.
//  The weights of a direction sum to one.
scalar conditionalWeight(const labelList& nodeIndex, const label dimi)
{
    const scalar p = 0.35 + 0.05*std::fmod(prefixKey(nodeIndex, dimi), 3.0);

    return nodeIndex[dimi] == 0 ? p : 1.0 - p;
}


//- Abscissa of the node of direction dimi, conditioned on the ones before
scalar conditionalAbscissa(const labelList& nodeIndex, const label dimi)
{
    const scalar centre = 0.2*dimi + 0.1*prefixKey(nodeIndex, dimi);
    const scalar spread = 1.0 + 0.1*dimi;

    return nodeIndex[dimi] == 0 ? centre - spread : centre + spread;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Invert a conditional quadrature one of whose abscissae is zero.
//
//  The distribution is a two by two conditional quadrature, so the method
//  represents it exactly and every moment has to come back. The first
//  direction carries its nodes at zero and at one, which is what a
//  distribution with a still fraction looks like, and a zero abscissa is an
//  ordinary node of a quadrature: it is the central node of any symmetric
//  distribution.
void testZeroAbscissa()
{
    Info<< "\n\nInverting a conditional quadrature with an abscissa at zero"
        << endl;

    // First direction: weights and abscissae, one of which is zero
    const scalarList w0({0.4, 0.6});
    const scalarList x0({0.0, 1.0});

    // Second direction, conditioned on the node of the first
    const List<scalarList> p1({{0.7, 0.3}, {0.25, 0.75}});
    const List<scalarList> y1({{-1.0, 2.0}, {0.5, 3.0}});

    labelListList momentOrders
    ({
        {0, 0}, {1, 0}, {2, 0}, {3, 0},
        {0, 1}, {1, 1},
        {0, 2}, {1, 2},
        {0, 3}, {1, 3}
    });

    const labelListList nodeIndexes({{0, 0}, {0, 1}, {1, 0}, {1, 1}});
    const labelList velocityIndexes({-1});

    multivariateMomentSet moments
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(2, supportType::R),
        SMALL,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scalar m = 0.0;

        forAll(w0, i)
        {
            scalar conditional = 0.0;

            forAll(p1[i], j)
            {
                conditional += p1[i][j]*pow(y1[i][j], momentOrder[1]);
            }

            m += w0[i]*pow(x0[i], momentOrder[0])*conditional;
        }

        moments(momentOrder) = m;
    }

    dictionary dict
    (
        IStringStream
        (
            "supports (\"R\" \"R\");"
            "basicQuadrature0 { univariateMomentInversion Gauss; }"
            "basicQuadrature1 { univariateMomentInversion Gauss; }"
        )()
    );

    multivariateMomentInversions::conditional inverter
    (
        dict, momentOrders, nodeIndexes, velocityIndexes
    );

    if (!inverter.invert(moments))
    {
        FatalErrorInFunction
            << "The inversion of the conditional quadrature failed." << nl
            << exit(FatalError);
    }

    const mappedScalarList& weights = inverter.weights();
    const mappedList<scalarList>& abscissae = inverter.abscissae();

    mappedList<scalar> reconstructed(momentOrders.size(), momentOrders, Zero);

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scalar m = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            m += weights(nodeIndex)
                *pow(abscissae(nodeIndex)[0], momentOrder[0])
                *pow(abscissae(nodeIndex)[1], momentOrder[1]);
        }

        reconstructed(momentOrder) = m;
    }

    checkMomentConservation
    (
        reconstructed,
        moments,
        momentOrders,
        1e-10,
        "conditional, with an abscissa at zero"
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "createFields.H"

    mappedList<scalarList> x
    (
        nNodes,
        nodeIndexes,
        scalarField(nDims, Zero)
    );
    mappedList<scalar> w(nNodes, nodeIndexes, 0.0);

    // A conditional quadrature of the shape the method builds, so that the
    // moments it is given are the moments of a distribution it represents
    forAll(nodeIndexes, nodei)
    {
        const labelList& nodeIndex = nodeIndexes[nodei];

        scalar weight = 1.0;

        forAll(nodeIndex, dimi)
        {
            weight *= conditionalWeight(nodeIndex, dimi);
            x(nodeIndex)[dimi] = conditionalAbscissa(nodeIndex, dimi);
        }

        w(nodeIndex) = weight;
    }

    Info<< "Original moments:" << endl;

    multivariateMomentSet moments
    (
        nMoments,
        momentOrders,
        List<supportType>(nDims, supportType::R),
        SMALL,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        moments(momentOrder) = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            scalar cmpt = w(nodeIndex);

            forAll(nodeIndex, dimi)
            {
                cmpt *= pow(x(nodeIndex)[dimi], momentOrder[dimi]);
            }

            moments(momentOrder) += cmpt;
        }

        Info<< "moment.";

        forAll(momentOrder, dimi)
        {
            Info<< momentOrder[dimi];
        }

        Info<< ": " << moments(momentOrder) << endl;
    }

    multivariateMomentInversions::conditional momentInverter
    (
        quadratureProperties, momentOrders, nodeIndexes, velocityIndexes
    );

    Info<< "\nInverting moments" << endl;

    momentInverter.invert(moments);

    Info<< "\nReconstructed moments:" << endl;

    const mappedScalarList& weights = momentInverter.weights();
    const mappedList<scalarList>& abscissae = momentInverter.abscissae();

    const mappedVectorList& velocityAbscissae =
        momentInverter.velocityAbscissae();

    mappedList<scalar> newMoments(nMoments, momentOrders);

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        newMoments(momentOrder) = 0.0;

        forAll(nodeIndexes, nodei)
        {
            const labelList& nodeIndex = nodeIndexes[nodei];

            scalar cmpt = weights(nodeIndex);
            label vi = 0;
            label si = 0;

            for(label dimi = 0; dimi < momentOrder.size(); dimi++)
            {
                if (vi < velocityIndexes.size() && velocityIndexes[vi] == dimi)
                {
                     cmpt *=
                        pow
                        (
                            velocityAbscissae(nodeIndex)[vi],
                            momentOrder[dimi]
                        );
                    vi++;
                }
                else
                {
                    cmpt *= pow(abscissae(nodeIndex)[si], momentOrder[dimi]);
                    si++;
                }
            }

            newMoments(momentOrder) += cmpt;
        }

    }

    // The distribution is a conditional quadrature of the shape the method
    // builds, so every moment of the set has to be reproduced
    checkMomentConservation
    (
        newMoments,
        moments,
        momentOrders,
        1e-10,
        "conditional, five dimensions"
    );


    testZeroAbscissa();

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
