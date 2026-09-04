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
    Test-TensorProductMomentInversion

Description
    Test the tensor product moment inversion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "mappedLists.H"
#include "supportType.H"
#include "TensorProductMomentInversion.H"
#include "Random.H"
#include "multivariateMomentTest.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Invert a moment set whose zero-order moment is below the smallest one
//  the inversion accepts.
//
//  Such a set carries no information, and the realizability check settles
//  it before any quadrature is built, so the inversion reports the failure
//  and leaves the quadrature at the zero it starts from. The caller is what
//  decides what to do with a cell like that, so it has to be told.
void testNegligibleMass
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
{
    Info<< "\n\nInverting a moment set of negligible mass" << endl;

    multivariateMomentSet moments
    (
        momentOrders.size(),
        momentOrders,
        List<supportType>(momentOrders[0].size(), supportType::R),
        SMALL,
        SMALL
    );

    forAll(momentOrders, mi)
    {
        moments(momentOrders[mi]) = 0.1*SMALL;
    }

    multivariateMomentInversions::TensorProduct inverter
    (
        dict, momentOrders, nodeIndexes, velocityIndexes
    );

    if (inverter.invert(moments))
    {
        FatalErrorInFunction
            << "The inversion of a moment set of negligible mass was"
            << " reported as a success." << nl
            << "    Zero-order moment: " << moments(momentOrders[0]) << nl
            << "    Smallest accepted: " << inverter.smallM0() << nl
            << exit(FatalError);
    }

    Info<< "  the inversion reported the failure" << endl;

    forAll(inverter.weights(), nodei)
    {
        if (mag(inverter.weights()[nodei]) > SMALL)
        {
            FatalErrorInFunction
                << "A failed inversion left a weight behind." << nl
                << "    Node: " << nodei << nl
                << "    Weight: " << inverter.weights()[nodei] << nl
                << exit(FatalError);
        }
    }

    Info<< "  the quadrature it leaves behind is null" << endl;
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
            x[nodei][dimi] = rndGen.sample01<scalar>();
        }
    }

    Info<< "Original moments:" << endl;

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

    multivariateMomentInversions::TensorProduct momentInverter
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

    // The tensor product quadrature has a node for every combination of the
    // univariate abscissae, and solves for its weights from the mixed
    // moments, so every moment of the set is reproduced
    checkMomentConservation
    (
        newMoments,
        moments,
        momentOrders,
        1e-10,
        "TensorProduct"
    );


    testNegligibleMass
    (
        quadratureProperties,
        momentOrders,
        nodeIndexes,
        velocityIndexes
    );

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
