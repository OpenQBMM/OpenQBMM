/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2023 Alberto Passalacqua
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
#include "TensorProductMomentInversion.H"
#include "Random.H"

using namespace Foam;

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

    forAll(x, nodei)
    {
        w[nodei] = scalar(rand())/scalar(RAND_MAX);

        forAll(x[nodei], dimi)
        {
            x[nodei][dimi] = scalar(rand())/scalar(RAND_MAX);
        }
    }

    Info<< "Original moments:" << endl;

    multivariateMomentSet moments(nMoments, momentOrders, "R", SMALL, SMALL);

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
                if (velocityIndexes[vi] == dimi)
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

        Info<< "moment.";

        forAll(momentOrder, dimi)
        {
            Info<< momentOrder[dimi];
        }

        Info<< ": " << newMoments(momentOrder)
            << ",\trel error: "
            << (mag(moments(momentOrder) 
                - newMoments(momentOrder))/moments(momentOrder))<< endl;
    }

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
