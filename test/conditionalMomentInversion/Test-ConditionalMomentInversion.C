/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
     \\/     M anipulation  |
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
    Test-ExtendedMomentInversion.C

Description
    Test the extendedMomentInversion class and its subclasses.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "nDimensionalMappedList.H"
#include "conditionalMomentInversion.H"
#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "createFields.H"

    label nMoments = 1;
    label nDims = 4;

    // Set number of nodes in each direction
    labelList nNodes(nDims);

    nNodes[0] = 1;
    nNodes[1] = 2;
    nNodes[2] = 3;
    nNodes[3] = 4;

    // Declaration of permutation ([0 1 2 3] is normal).
    labelList permutation(nDims);

    permutation[0] = 0;
    permutation[1] = 1;
    permutation[2] = 2;
    permutation[3] = 3;


    labelList nNodesP = nNodes;

    for (label dimi = 0; dimi < nDims; dimi++)
    {
        nNodesP[dimi] = nNodes[permutation[dimi]];
    }

    // Set number of moments corresponding to 2*nNodes[i] in each direction
    // (Worst case senario)
    for (int i = 0; i < nDims; i++)
    {
        nMoments *= 2*nNodes[i];
    }

    List<word> support(nDims, "RPlus");

    PtrList<nDimensionalMappedList<scalar>> x(nDims);
    PtrList<nDimensionalMappedList<scalar>> w(nDims);

    forAll(x, i)
    {
        x.set
        (
            i,
            new nDimensionalMappedList<scalar>
            (
                i + 1,
                nNodes
            )
        );

        w.set
        (
            i,
            new nDimensionalMappedList<scalar>
            (
                i + 1,
                nNodes
            )
        );

        forAll(x[i], j)
        {
            x[i].set
            (
                j,
                new scalar(0.0)
            );

            w[i].set
            (
                j,
                new scalar(0.0)
            );
        }
    }

    forAll(x, dimi)
    {
        forAll(x[dimi], cmpti)
        {
            x[dimi][cmpti] = scalar(rand())/scalar(RAND_MAX);
            w[dimi][cmpti] = scalar(rand())/scalar(RAND_MAX);
        }
    }

    Map<label> map(nMoments);
    Map<label> mapP(nMoments);

    labelList pos(nDims);
    labelList posP(nDims);

    label mi = 0;

    for (label i = 0; i < 2*nNodes[0]; i++)
    {
        for (label j = 0; j < 2*nNodes[1]; j++)
        {
            for (label k = 0; k < 2*nNodes[2]; k++)
            {
                for (label l = 0; l < 2*nNodes[3]; l++)
                {
                    pos[0] = i;
                    pos[1] = j;
                    pos[2] = k;
                    pos[3] = l;

                    map.insert
                    (
                        Foam::nDimensionalMappedList<scalar>::listToLabel
                        (
                            pos
                        ),
                        mi
                    );

                    posP[0] = pos[permutation[0]];
                    posP[1] = pos[permutation[1]];
                    posP[2] = pos[permutation[2]];
                    posP[3] = pos[permutation[3]];

                    mapP.insert
                    (
                        Foam::nDimensionalMappedList<scalar>::listToLabel
                        (
                            posP
                        ),
                        mi
                    );

                    mi++;
                }
            }
        }
    }

    Info<< "Original moments:" << endl;

    nDimensionalMappedList<scalar> moments(nMoments, nDims, map);
    nDimensionalMappedList<scalar> momentsP(nMoments, nDims, mapP);

    forAll(moments, mi)
    {
        moments.set
        (
            mi,
            new scalar(0.0)
        );

        momentsP.set
        (
            mi,
            new scalar(0.0)
        );
    }

    scalar sum = 0.0;

    for (label l = 0; l < 2*nNodes[0]; l++)
    {
        pos[0] = l;

        for (label m = 0; m < 2*nNodes[1]; m++)
        {
            pos[1] = m;

            for (label n = 0; n < 2*nNodes[2]; n++)
            {
                pos[2] = n;

                for (label nn = 0; nn < 2*nNodes[3]; nn++)
                {
                    pos[3] = nn;

                    sum = 0.0;

                    for (label i = 0; i < nNodes[0]; i++)
                    {
                        for (label j = 0; j < nNodes[1]; j++)
                        {
                            for (label k = 0; k < nNodes[2]; k++)
                            {
                                for (label kk = 0; kk < nNodes[3]; kk++)
                                {
                                    sum += w[0](i)*w[1](i, j)*w[2](i, j, k)
                                        *w[3](i, j, k, kk)
                                        *pow(x[0](i), l)*pow(x[1](i, j), m)
                                        *pow(x[2](i, j, k), n)
                                        *pow(x[3](i, j, k, kk), nn);
                                }
                            }
                        }
                    }

                    moments(l, m, n, nn) = sum;

                    Info<< "moment." << l << m << n << nn <<": "
                        << moments(l, m, n, nn) << endl;
                }
            }
        }
    }

    // Premutation 213
    for (label l = 0; l < 2*nNodes[0]; l++)
    {
        for (label m = 0; m < 2*nNodes[1]; m++)
        {
            for (label n = 0; n < 2*nNodes[2]; n++)
            {
                for (label nn = 0; nn < 2*nNodes[3]; nn++)
                {
                    pos[0] = l;
                    pos[1] = m;
                    pos[2] = n;
                    pos[3] = nn;

                    posP[0] = pos[permutation[0]];
                    posP[1] = pos[permutation[1]];
                    posP[2] = pos[permutation[2]];
                    posP[3] = pos[permutation[3]];

                    momentsP(posP) = moments(pos);
                }
            }
        }
    }

    conditionalMomentInversion momentInverter
    (
        quadratureProperties, nMoments, mapP, nNodesP, support
    );

    Info<< "\nInverting moments" << endl;

    momentInverter.invert(momentsP);

    Info<< "\nReconstructed moments:" << endl;

    const PtrList<nDimensionalMappedList<scalar> >& weights =
        momentInverter.weights();

    const PtrList<nDimensionalMappedList<scalar> >& abscissae =
        momentInverter.abscissae();

    for (label l = 0; l < 2*nNodesP[0]; l++)
    {
        pos[0] = l;

        for (label m = 0; m < 2*nNodesP[1]; m++)
        {
            pos[1] = m;

            for (label n = 0; n < 2*nNodesP[2]; n++)
            {
                pos[2] = n;

                for (label nn = 0; nn < 2*nNodesP[3]; nn++)
                {
                    pos[3] = nn;
                    sum = 0.0;

                    for (label i = 0; i < nNodesP[0]; i++)
                    {
                        for (label j = 0; j < nNodesP[1]; j++)
                        {
                            for (label k = 0; k < nNodesP[2]; k++)
                            {
                                for (label kk = 0; kk < nNodesP[3]; kk++)
                                {
                                    sum += weights[0](i)*weights[1](i, j)
                                        *weights[2](i, j, k)
                                        *weights[3](i, j, k, kk)
                                        *pow(abscissae[0](i), l)
                                        *pow(abscissae[1](i, j), m)
                                        *pow(abscissae[2](i, j, k), n)
                                        *pow(abscissae[3](i, j, k, kk), nn);
                                }
                            }
                        }
                    }

                    momentsP(l, m, n, nn) = sum;
                }
            }
        }
    }

    for (label l = 0; l < 2*nNodes[0]; l++)
    {
        for (label m = 0; m < 2*nNodes[1]; m++)
        {
            for (label n = 0; n < 2*nNodes[2]; n++)
            {
                for (label nn = 0; nn < 2*nNodes[3]; nn++)
                {
                    pos[0] = l;
                    pos[1] = m;
                    pos[2] = n;
                    pos[3] = nn;

                    posP[0] = pos[permutation[0]];
                    posP[1] = pos[permutation[1]];
                    posP[2] = pos[permutation[2]];
                    posP[3] = pos[permutation[3]];

                    moments(pos) = momentsP(posP);

                    Info<< "moment." << l << m << n << nn << ": "
                        << moments(pos) << endl;
                }
            }
        }
    }

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
