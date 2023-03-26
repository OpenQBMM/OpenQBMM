/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
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
#include "scalarList.H"
#include "scalarMatrices.H"
#include "univariateMomentSet.H"
#include "extendedMomentInversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<< "Reading quadratureProperties\n" << endl;
    dictionary quadratureProperties(IFstream("quadratureProperties")());

    Info << "Testing extendedMomentInversion\n" << endl;

    label nMoments = 5;
    word support = "RPlus";
    univariateMomentSet moments(nMoments, support, SMALL, SMALL);

    // Dirac delta function
//     moments[0] = 1.0;
//     moments[1] = 1.0;
//     moments[2] = 1.0;
//     moments[3] = 1.0;
//     moments[4] = 1.0;

//  Valid moment set
   moments[0] = 1.0;
   moments[1] = 2.708217669;
   moments[2] = 8.951330468;
   moments[3] = 35.95258119;
   moments[4] = 174.4370267;

//  Unrealizable moment star
//     moments[0] = 0.567128698550116;
//     moments[1] = 0.659798044636756;
//     moments[2] = 0.796168501018439;
//     moments[3] = 1.0;
//     moments[4] = 1.3103698092;

//  Set of moment with sigma = 0 as root
//     moments[0] = 1.0;
//     moments[1] = 1.6487212707;
//     moments[2] = 2.7182818285;
//     moments[3] = 4.4816890703;
//     moments[4] = 7.3890560989;

// Set of moments with multiple roots not bracketed by sigmaMax
//     moments[0] = 1.0;
//     moments[1] = 1.055795897;
//     moments[2] = 1.247672408;
//     moments[3] = 1.670093417;
//     moments[4] = 2.578636894;
//     moments[5] = 4.688985807;
//     moments[6] = 10.22526492;
//     moments[7] = 27.0370458224;
//     moments[8] = 86.9534420717;

// Moment set on edge of moment space
//     moments[0] = 3.125e12;
//     moments[1] = 6.25e6;
//     moments[2] = 12.5;
//     moments[3] = 2.5e-5;
//     moments[4] = 5.0e-11;

//     moments[0] = 1.0;
//     moments[1] = 2.0;
//     moments[2] = 4.0;
//     moments[3] = 8.0;
//     moments[4] = 16.0;

//     moments[0] = 0.9996;
//     moments[1] =  0.99970396842;
//     moments[2] = 0.999834960421;
//     moments[3] = 1;
//     moments[4] = 1.00020793684;


// Moment set provided by Frederique Laurent-Negre
// M_k = 1/k; k = 1, 2*N + 1
// Test with 11 moments - All should be reproduced
//        for (label mI = 1; mI < nMoments + 1; mI++)
//        {
//            moments[mI - 1] = 1.0/scalar(mI);
//        }
//
//        moments[9] = 0;

// Moment set provided by Frederique Laurent-Negre
// Last moment not preserved
//    moments[0] = 1.0;
//    moments[1] = 1.0;
//    moments[2] = 1.1;
//    moments[3] = 1.41;
//    moments[4] = 2.371;
//    moments[5] = 5.4501;
//    moments[6] = 15.75531;

    Info << setprecision(16);
    Info << "Input moments\n" << endl;

    for (label momentI = 0; momentI < nMoments; momentI++)
    {
        Info << "Moment " << momentI << " = " << moments[momentI] << endl;
    }

    Info << endl;

    autoPtr<extendedMomentInversion> EQMOM
    (
        extendedMomentInversion::New
        (
            quadratureProperties,
            nMoments,
            readLabel(quadratureProperties.lookup("nSecondaryNodes"))
        )
    );

    Info << "\nInverting moments.\n" << endl;

    EQMOM->invert(moments);

    Info << "Sigma = " << EQMOM->sigma() << endl;
    Info << "\nExtracting secondary quadrature." << endl;
    Info << "\nRecovering secondary weights and abscissae." << endl;

    const scalarList& pWeights(EQMOM->primaryWeights());
    const scalarList& pAbscissae(EQMOM->primaryAbscissae());
    const scalarRectangularMatrix& sWeights(EQMOM->secondaryWeights());
    const scalarRectangularMatrix& sAbscissae(EQMOM->secondaryAbscissae());

    Info << "\nStoring quadrature." << endl;

    OFstream outputFile("./secondaryQuadrature");

    label nPrimaryNodes = EQMOM->nPrimaryNodes();
    label nSecondaryNodes = EQMOM->nSecondaryNodes();

    for (label pNodeI = 0; pNodeI < nPrimaryNodes; pNodeI++)
    {
        outputFile << "Primary node " << pNodeI
            << "\nPrimary weight = " << pWeights[pNodeI]
            << "\nPrimary abscissa = " << pAbscissae[pNodeI] << endl;

        outputFile << "\nSecondary nodes" << endl;

        for (label sNodeI = 0; sNodeI < nSecondaryNodes; sNodeI++)
        {
            outputFile << sWeights[pNodeI][sNodeI] << ", "
                << sAbscissae[pNodeI][sNodeI] << endl;
        }

        outputFile << "\n\n";
    }

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
