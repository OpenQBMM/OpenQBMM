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

    label nMoments = 10;
    label nNodes = 9;

    labelListList indicies(nNodes, labelList(2,0));
    indicies[0] = {1,1};
    indicies[1] = {2,1};
    indicies[2] = {1,2};
    indicies[3] = {2,2};
    indicies[4] = {3,1};
    indicies[5] = {1,3};
    indicies[6] = {2,3};
    indicies[7] = {3,2};
    indicies[8] = {3,3};

    mappedList<scalar> w(nNodes, indicies);
    mappedList<vector> u(nNodes, indicies);

//     for(label i = 1; i <= 3; i++)
//     {
//         for(label j = 1; j <= 3; j++)
//         {
//             w(i,j) = scalar(rand())/scalar(RAND_MAX);
//             u(i,j) =
//                 vector
//                 (
//                     0.5,//scalar(rand())/scalar(RAND_MAX),
//                     0.5,//scalar(rand())/scalar(RAND_MAX),
//                     0
//                 )*10 - vector(5.0, 5.0, 0);
//         }
//     }
    w(1,1) = 0; u(1,1) = vector(0,0,0);
    w(1,2) = 0; u(1,2) = vector(0,0,0);
    w(1,3) = 0; u(1,3) = vector(0,0,0);
    w(2,1) = 0; u(2,1) = vector(0,0,0);
    w(2,2) = 0; u(2,2) = vector(0,0,0);
    w(2,3) = 0; u(2,3) = vector(0,0,0);
    w(3,1) = 0; u(3,1) = vector(0,0,0);
    w(3,2) = 0; u(3,2) = vector(0,0,0);
    w(3,3) = 0; u(3,3) = vector(0,0,0);


    labelListList mIndicies(10, labelList(2,0));
    mIndicies[0] = {0,0};
    mIndicies[1] = {1,0};
    mIndicies[2] = {0,1};
    mIndicies[3] = {2,0};
    mIndicies[4] = {0,2};
    mIndicies[5] = {3,0};
    mIndicies[6] = {0,3};
    mIndicies[7] = {4,0};
    mIndicies[8] = {0,4};
    mIndicies[9] = {1,1};

    Info<< "Original moments:" << endl;

    multivariateMomentSet moments(nMoments, mIndicies, "R");

//     for (label i = 1; i <= 3; i++)
//     {
//         for (label j = 1; j <= 3; j++)
//         {
//             forAll(moments, mi)
//             {
//                 moments[mi] +=
//                     w(i,j)
//                    *pow(u(i,j).x(), mIndicies[mi][0])
//                    *pow(u(i,j).y(), mIndicies[mi][1]);
//             }
//         }
//     }
//     moments(0,0) = 5.43408e-15;
//     moments(0,1) = 2.72019e-15;
//     moments(0,2) = 1.3627e-15;
//     moments(0,3) = 6.83773e-16;
//     moments(0,4) = 3.43956e-16;
//     moments(1,0) = -5.44038e-15;
//     moments(1,1) = -2.7254e-15;
//     moments(2,0) = 5.45079e-15;
//     moments(3,0) = -5.47081e-15;
//     moments(4,0) = 5.5033e-15;
    moments(0,0) = 0.000182424;
    moments(1,0) = 0.000181355;
    moments(0,1) = 0.000182424;
    moments(2,0) = 0.000182362;
    moments(1,1) = 0.000181355;
    moments(0,2) = 0.000182424;
    moments(3,0) = 0.000181189;
    moments(0,3) = 0.000182424;
    moments(4,0) = 0.000182213;
    moments(0,4) = 0.000182424;

    forAll(mIndicies, mi)
    {
        Info<< "moment." << mIndicies[mi][0] << mIndicies[mi][1]<<": "
            << moments[mi] <<endl;
    }

    hyperbolicConditionalMomentInversion momentInverter
    (
        quadratureProperties, 2
    );

    Info<< "\nInverting moments" << endl;

    momentInverter.invert(moments);

    Info<< "\nReconstructed moments:" << endl;

    const mappedList<scalar>& weights = momentInverter.weights();
    const mappedList<vector>& abscissae = momentInverter.abscissae();
    multivariateMomentSet momentsNew(nMoments, mIndicies, "R");


    for (label i = 1; i <= 3; i++)
    {
        for (label j = 1; j <= 3; j++)
        {
            forAll(mIndicies, mi)
            {
                momentsNew[mi] +=
                    weights(i,j)
                   *pow(abscissae(i,j).x(), mIndicies[mi][0])
                   *pow(abscissae(i,j).y(), mIndicies[mi][1]);
            }
            Info<<i<<j<<"-u: " << abscissae(i,j) << nl
                    <<"w: " <<weights(i,j)<<endl;
        }
    }

    forAll(moments, mi)
    {
        Info<< "moment." << mIndicies[mi][0] << mIndicies[mi][1]<<": "
            << momentsNew[mi] <<endl;
    }
    Info << nl << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
