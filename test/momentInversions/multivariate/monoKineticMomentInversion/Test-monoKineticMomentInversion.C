/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2026 Alberto Passalacqua
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
    Test-monokineticMomentInversion

Description
    Test the monoKinetic multivariate moment inversion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "mappedLists.H"
#include "supportType.H"
#include "monoKineticMomentInversion.H"
#include "Random.H"
#include "multivariateMomentTest.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Invert a monokinetic distribution whose dimensions are ordered velocity
//  first and size last.
//
//  A monokinetic distribution carries one velocity per size node, so it is
//  its own quadrature and the inversion has to reproduce every moment of
//  it. The point of the permuted order is the size: an inversion that
//  reads velocityIndexes to find out which dimension is which handles it,
//  while one that assumes the size is the first dimension reads the moments
//  of the wrong orders.
void testPermutedDimensions(const dictionary& dict)
{
    Info<< "\n\nInverting a distribution with the size as the last"
        << " dimension" << endl;

    const label nSizeNodes = 3;
    const label nSizeMoments = 2*nSizeNodes;

    // One velocity per size node
    const scalarList ws({0.5, 0.75, 1.0});
    const scalarList xs({0.4, 1.1, 1.8});
    const scalarList us({-0.6, 0.3, 1.2});
    const scalarList vs({0.9, -0.4, 0.2});

    // The dimensions are (u, v, size)
    labelListList momentOrders;

    for (label k = 0; k < nSizeMoments; k++)
    {
        momentOrders.append(labelList({0, 0, k}));
    }

    for (label k = 0; k < nSizeNodes; k++)
    {
        momentOrders.append(labelList({1, 0, k}));
        momentOrders.append(labelList({0, 1, k}));
    }

    const labelListList nodeIndexes({{0}, {1}, {2}});
    const labelList velocityIndexes({0, 1});

    multivariateMomentSet moments
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(3, supportType::R),
        SMALL,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scalar m = 0.0;

        forAll(ws, nodei)
        {
            m += ws[nodei]
                *pow(us[nodei], momentOrder[0])
                *pow(vs[nodei], momentOrder[1])
                *pow(xs[nodei], momentOrder[2]);
        }

        moments(momentOrder) = m;
    }

    multivariateMomentInversions::monoKinetic inverter
    (
        dict, momentOrders, nodeIndexes, velocityIndexes
    );

    if (!inverter.invert(moments))
    {
        FatalErrorInFunction
            << "The inversion of the permuted distribution failed." << nl
            << exit(FatalError);
    }

    const mappedScalarList& weights = inverter.weights();
    const mappedList<scalarList>& sizeAbscissae = inverter.abscissae();
    const mappedVectorList& velocityAbscissae = inverter.velocityAbscissae();

    mappedList<scalar> reconstructed(momentOrders.size(), momentOrders, Zero);

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];

        scalar m = 0.0;

        forAll(nodeIndexes, nodei)
        {
            m += weights[nodei]
                *pow(velocityAbscissae[nodei][0], momentOrder[0])
                *pow(velocityAbscissae[nodei][1], momentOrder[1])
                *pow(sizeAbscissae[nodei][0], momentOrder[2]);
        }

        reconstructed(momentOrder) = m;
    }

    checkMomentConservation
    (
        reconstructed,
        moments,
        momentOrders,
        1e-10,
        "monoKinetic, size as the last dimension"
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

    // A fixed seed, so that the moments the inversion is asked for are the
    // same on every run and on every machine
    Random rndGen(20260904);

    forAll(x, nodei)
    {
        w[nodei] = rndGen.sample01<scalar>();

        forAll(x[nodei], dimi)
        {
            if (dimi == 0)
            {
                x[nodei][dimi] = rndGen.sample01<scalar>();
            }
            else
            {
                x[nodei][dimi] = -0.5 + rndGen.sample01<scalar>();
            }
        }
    }

    Info<< "Initial abscissae: \n" << x << endl;
    Info<< "Initial weights: \n" << w << endl;

    Info<< "\nOriginal moments:" << endl;

    multivariateMomentSet moments
    (
        nMoments,
        momentOrders,
        List<supportType>(momentOrders[0].size(), supportType::R),
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

            cmpt *= pow(x[nodei][0], momentOrder[0]);

            forAll(velocityIndexes, dimi)
            {
                cmpt *=
                    pow
                    (
                        x[nodei][dimi + 1],
                        momentOrder[dimi + 1]
                    );
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

    multivariateMomentInversions::monoKinetic momentInverter
    (
        quadratureProperties, momentOrders, nodeIndexes, velocityIndexes
    );

    Info<< "\nInverting moments" << endl;

    momentInverter.invert(moments);

    const mappedScalarList& weights = momentInverter.weights();
    const mappedList<scalarList>& sizeAbscissae = momentInverter.abscissae();

    const mappedVectorList& velocityAbscissae =
        momentInverter.velocityAbscissae();

    Info << "Weights: " << weights << endl;
    Info << "Size abscissae: " << sizeAbscissae << endl;
    Info << "Velocity abscissae: " << velocityAbscissae << endl;

    Info<< "\nReconstructed moments:" << endl;

    mappedList<scalar> newMoments(nMoments, momentOrders);

    forAll(momentOrders, mi)
    {
        const labelList& momentOrder = momentOrders[mi];
        newMoments(momentOrder) = 0.0;

        forAll(nodeIndexes, nodei)
        {
            scalar cmpt = weights[nodei];
            cmpt *= pow(sizeAbscissae[nodei][0], momentOrder[0]);

            forAll(velocityIndexes, dimi)
            {
                cmpt *=
                    pow
                    (
                        velocityAbscissae[nodei][dimi],
                        momentOrder[dimi + 1]
                    );
            }

            newMoments(momentOrder) += cmpt;
        }

    }

    // monoKinetic assigns one velocity to each size node, so the quadrature
    // it builds reproduces every moment of the set it is given
    checkMomentConservation
    (
        newMoments,
        moments,
        momentOrders,
        1e-10,
        "monoKinetic"
    );

    testPermutedDimensions(quadratureProperties);

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
