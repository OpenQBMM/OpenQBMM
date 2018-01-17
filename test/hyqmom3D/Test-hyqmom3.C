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
#include "mappedList.H"
#include "hyperbolicConditionalMomentInversion.H"
#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "createFields.H"

    label nMoments = 16;

    labelListList indicies(27, labelList(3,0));
    indicies[0] = {1, 1, 1};
    indicies[1] = {1, 1, 2};
    indicies[2] = {1, 1, 3};
    indicies[3] = {1, 2, 1};
    indicies[4] = {1, 2, 2};
    indicies[5] = {1, 2, 3};
    indicies[6] = {1, 3, 1};
    indicies[7] = {1, 3, 2};
    indicies[8] = {1, 3, 3};
    indicies[9] = {2, 1, 1};
    indicies[10] = {2, 1, 2};
    indicies[11] = {2, 1, 3};
    indicies[12] = {2, 2, 1};
    indicies[13] = {2, 2, 2};
    indicies[14] = {2, 2, 3};
    indicies[15] = {2, 3, 1};
    indicies[16] = {2, 3, 2};
    indicies[17] = {2, 3, 3};
    indicies[18] = {3, 1, 1};
    indicies[19] = {3, 1, 2};
    indicies[20] = {3, 1, 3};
    indicies[21] = {3, 2, 1};
    indicies[22] = {3, 2, 2};
    indicies[23] = {3, 2, 3};
    indicies[24] = {3, 3, 1};
    indicies[25] = {3, 3, 2};
    indicies[26] = {3, 3, 3};

    mappedList<scalar> w(27, indicies);
    mappedList<vector> u(27, indicies);

    for (label i = 1; i <= 3; i++)
    {
        for(label j = 1; j <= 3; j++)
        {
            for(label k = 1; k <= 3; k++)
            {
//         label i = indicies[ind][0];
//         label j = indicies[ind][1];
//         label k = indicies[ind][2];

                w(i,j,k) = scalar(rand())/scalar(RAND_MAX);
                u(i,j,k) =
                    vector
                    (
                        scalar(rand())/scalar(RAND_MAX),
                        scalar(rand())/scalar(RAND_MAX),
                        scalar(rand())/scalar(RAND_MAX)
                    )*10 - vector(5.0, 5.0, 5.0);

                Info<<i<<j<<k<<": " << w(i,j,k)
                    <<"\t"<<u(i,j,k)<<endl;
            }
        }
    }

    labelListList mIndicies(16, labelList(3,0));
    mIndicies[0] = {0, 0, 0};
    mIndicies[1] = {1, 0, 0};
    mIndicies[2] = {0, 1, 0};
    mIndicies[3] = {0, 0, 1};
    mIndicies[4] = {2, 0, 0};
    mIndicies[5] = {1, 1, 0};
    mIndicies[6] = {1, 0, 1};
    mIndicies[7] = {0, 2, 0};
    mIndicies[8] = {0, 1, 1};
    mIndicies[9] = {0, 0, 2};
    mIndicies[10] = {3, 0, 0};
    mIndicies[11] = {0, 3, 0};
    mIndicies[12] = {0, 0, 3};
    mIndicies[13] = {4, 0, 0};
    mIndicies[14] = {0, 4, 0};
    mIndicies[15] = {0, 0, 4};

    Info<< "Original moments:" << endl;

    multivariateMomentSet moments(nMoments, mIndicies, "R");
    multivariateMomentSet momentsNew(nMoments, mIndicies, "R");

    for (label i = 1; i <= 3; i++)
    {
        for (label j = 1; j <= 3; j++)
        {
            for (label k = 1; k <= 3; k++)
            {
                forAll(moments, mi)
                {
                    moments[mi] +=
                        w(i,j,k)
                       *pow(u(i,j,k).x(), mIndicies[mi][0])
                       *pow(u(i,j,k).y(), mIndicies[mi][1])
                       *pow(u(i,j,k).z(), mIndicies[mi][2]);
                }
            }
        }
    }
    forAll(mIndicies, mi)
    {
        Info<< "moment."
            << mIndicies[mi][0]
            << mIndicies[mi][1]
            << mIndicies[mi][2]
            << ": " << moments[mi] <<endl;
    }

    hyperbolicConditionalMomentInversion momentInverter
    (
        quadratureProperties, 3
    );

    Info<< "\nInverting moments" << endl;
    momentInverter.invert(moments);

    Info<< "\nReconstructed moments:" << endl;

    const mappedList<scalar>& weights = momentInverter.weights();
    const mappedList<vector>& abscissae = momentInverter.abscissae();

    for (label i = 1; i <= 3; i++)
    {
        for (label j = 1; j <= 3; j++)
        {
            for (label k = 1; k <= 3; k++)
            {
                forAll(moments, mi)
                {
                    momentsNew[mi] +=
                        weights(i,j,k)
                       *pow(abscissae(i,j,k).x(), mIndicies[mi][0])
                       *pow(abscissae(i,j,k).y(), mIndicies[mi][1])
                       *pow(abscissae(i,j,k).z(), mIndicies[mi][2]);
                }
//                 Info<<i<<j<<k<<": " << weights(i,j,k)
//                     <<"\t"<<abscissae(i,j,k)<<endl;
            }
        }
    }

    forAll(moments, mi)
    {
        Info<< "moment."
            << mIndicies[mi][0]
            << mIndicies[mi][1]
            << mIndicies[mi][2]
            << ": " << momentsNew[mi] <<endl;
    }
    Info << nl << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
