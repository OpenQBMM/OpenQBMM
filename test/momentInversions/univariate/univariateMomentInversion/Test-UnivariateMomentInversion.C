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
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    Test-UnivariateMomentInversion

Description
    Test univariateMomentInversion class and methods.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "IFstream.H"
#include "OFstream.H"
#include "scalarMatrices.H"
#include "IOdictionary.H"
#include "univariateMomentSet.H"
#include "univariateMomentInversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void compareQuadrature
(
    scalarList& expectedWeights, 
    scalarList& computedWeights,
    scalarList& expectedAbscissae,
    scalarList& computedAbscissae
)
{
    Info<< "\n" << endl;

    if (expectedWeights.size() != computedWeights.size())
    {
        FatalErrorInFunction
            << "Weights vectors have different size: " 
            << endl
            << "  Size of computed weights vector: " << computedWeights.size()
            << endl
            << "  Size of expected weights vector: " << expectedWeights.size()
            << endl
            << exit(FatalError);
    }

    if (expectedAbscissae.size() != computedAbscissae.size())
    {
        FatalErrorInFunction
            << "Abscissae vectors have different size: " 
            << endl
            << "  Size of computed abscissae vector: " 
            << computedAbscissae.size() << endl
            << "  Size of expected abscissae vector: " 
            << expectedAbscissae.size()
            << endl
            << exit(FatalError);
    }

    forAll(expectedWeights, weightsi)
    {
        scalar magDiff 
            = mag(expectedWeights[weightsi] - computedWeights[weightsi]);

        Info<< "  expectedWeights[" << weightsi << "] = " 
            << expectedWeights[weightsi] 
            << ", computedWeights[" << weightsi << "] = " 
            << computedWeights[weightsi] 
            << endl;

        if (magDiff >= SMALL)
        {
            FatalErrorInFunction
                << "Values of the quadrature weights do not match: " 
                << endl
                << "  Position: " << weightsi
                << endl
                << "  Expected weight value " << expectedAbscissae[weightsi]
                << endl
                << "  Computed weight value " << computedAbscissae[weightsi]
                << endl
                << "  Magnitude of the difference " << magDiff
                << endl
                << exit(FatalError);
        }
    }

    Info<< "\nQuadrature weights match.\n" << endl;

    forAll(expectedAbscissae, abscissai)
    {
        scalar magDiff 
            = mag(expectedAbscissae[abscissai] - computedAbscissae[abscissai]);

        Info<< "  expectedAbscissae[" << abscissai << "] = " 
            << expectedAbscissae[abscissai] 
            << ", computedAbscissae[" << abscissai << "] = " 
            << computedAbscissae[abscissai] 
            << endl;

        if (magDiff >= SMALL)
        {
            FatalErrorInFunction
                << "Values of the quadrature abscissae do not match: " 
                << endl
                << "  Position: " << abscissai
                << endl
                << "  Expected abscissa value " << expectedAbscissae[abscissai]
                << endl
                << "  Computed abscissa value " << computedAbscissae[abscissai]
                << endl
                << "  Magnitude of the difference " << magDiff
                << endl
                << exit(FatalError);
        }
    }

    Info<< "\nQuadrature abscissae match.\n" << endl;
}

void showInputMoments(univariateMomentSet& inputMoments, string quadratureName)
{
    Info<< "\nTesting " << quadratureName << " quadrature\n" << endl;

    forAll(inputMoments, mi)
    {
        Info<< "  inputMoments[" << mi << "] = " << inputMoments[mi] << endl;
    }

    Info<< "\n" << endl;
}

void testQuadrature
(
    univariateMomentSet& m, 
    dictionary& dict,
    scalarList& expectedWeights,
    scalarList& expectedAbscissae,
    string quadratureName
)
{
    showInputMoments(m, quadratureName);

    autoPtr<univariateMomentInversion> inversion
    (
        univariateMomentInversion::New(dict)
    );

    inversion().invert(m, 0, 1);

    scalarList weights(inversion().weights());
    scalarList abscissae(inversion().abscissae());

    compareQuadrature(expectedWeights, weights, expectedAbscissae, abscissae);

    m.update(weights, abscissae);

    Info<< "\nMoments computed from quadrature\n" << endl;

    forAll(m, mi)
    {
        Info<< "  Moment " << mi << " = " << m[mi] << endl;
    }
}

int main(int argc, char *argv[])
{
    Info<< setprecision(16);

    Info<< "Testing univariateMomentInversion\n" << endl;
    Info<< "---------------------------------\n" << endl;

    Info<< "Reading quadraturePropertiesGauss\n" << endl;

    dictionary quadraturePropertiesGauss
    (
        IFstream("quadraturePropertiesGauss")()
    );

    Info<< "Reading quadraturePropertiesRadau\n" << endl;

    dictionary quadraturePropertiesRadau
    (
        IFstream("quadraturePropertiesRadau")()
    );

    Info<< "Reading quadraturePropertiesLobatto\n" << endl;

    dictionary quadraturePropertiesLobatto
    (
        IFstream("quadraturePropertiesLobatto")()
    );

    // Test 2 - m = (1, 1, 1, 1)
    scalarList inputMoments1(4, 1.0);
    univariateMomentSet mGaussTest1(inputMoments1, "RPlus");
    scalarList expectedWeightsTest1(1, 1.0);
    scalarList expectedAbscissaeTest1(1, 1.0);

    testQuadrature
    (
        mGaussTest1,
        quadraturePropertiesGauss,
        expectedWeightsTest1,
        expectedAbscissaeTest1,
        "Gauss"
    );

    // Test 2 - m_{i-1} = 1/m_i, with i = 1, ..., nMoments - 1
    scalarList inputMoments2(10);

    for (label mi = 1; mi < inputMoments2.size() + 1; mi++)
    {
        inputMoments2[mi - 1] = 1.0/scalar(mi);
    }

    univariateMomentSet mGaussTest2(inputMoments2, "RPlus");
    
    scalarList expectedWeightsTest2(5);

    expectedWeightsTest2[0] = 0.1184634425280107;
    expectedWeightsTest2[1] = 0.2393143352499202;
    expectedWeightsTest2[2] = 0.2844444444446074;
    expectedWeightsTest2[3] = 0.2393143352495329;
    expectedWeightsTest2[4] = 0.1184634425279286;

    scalarList expectedAbscissaeTest2(5);

    expectedAbscissaeTest2[0] = 0.04691007703061419;
    expectedAbscissaeTest2[1] = 0.2307653449471704;
    expectedAbscissaeTest2[2] = 0.5000000000002638;
    expectedAbscissaeTest2[3] = 0.7692346550530993;
    expectedAbscissaeTest2[4] = 0.9530899229694018;

    testQuadrature
    (
        mGaussTest2,
        quadraturePropertiesGauss,
        expectedWeightsTest2,
        expectedAbscissaeTest2,
        "Gauss"
    );

    // Test 3 - m_{i-1} = 1/m_i, with i = 1, ..., nMoments - 1, Gauss-Lobatto
    scalarList inputMoments3(11);

    for (label mi = 1; mi < inputMoments3.size() + 1; mi++)
    {
        inputMoments3[mi - 1] = 1.0/scalar(mi);
    }

    univariateMomentSet mGaussTest3(inputMoments3, "RPlus", 1);
    
    scalarList expectedWeightsTest3(6);

    expectedWeightsTest3[0] = 0.02777777777722189;
    expectedWeightsTest3[1] = 0.1598203766073751;
    expectedWeightsTest3[2] = 0.2426935942324191;
    expectedWeightsTest3[3] = 0.2604633915958309;
    expectedWeightsTest3[4] = 0.2084506671586522;
    expectedWeightsTest3[5] = 0.1007941926285008;

    scalarList expectedAbscissaeTest3(6);

    expectedAbscissaeTest3[0] = 0.0;
    expectedAbscissaeTest3[1] = 0.09853508579692917;
    expectedAbscissaeTest3[2] = 0.30453572664167;
    expectedAbscissaeTest3[3] = 0.562025189747392;
    expectedAbscissaeTest3[4] = 0.8019865821232577;
    expectedAbscissaeTest3[5] = 0.9601901429478156;

    testQuadrature
    (
        mGaussTest3,
        quadraturePropertiesRadau,
        expectedWeightsTest3,
        expectedAbscissaeTest3,
        "Gauss-Radau"
    );

    // Test 4 - m_{i-1} = 1/m_i, with i = 1, ..., nMoments - 1, Gauss-Radau
    scalarList inputMoments4(12);

    for (label mi = 1; mi < inputMoments4.size() + 1; mi++)
    {
        inputMoments4[mi - 1] = 1.0/scalar(mi);
    }

    univariateMomentSet mGaussTest4(inputMoments4, "RPlus", 2);
    
    scalarList expectedWeightsTest4(7);

    expectedWeightsTest4[0] = 0.02380952380528316;
    expectedWeightsTest4[1] = 0.1384130236614941;
    expectedWeightsTest4[2] = 0.215872690592007;
    expectedWeightsTest4[3] = 0.2438095238141702;
    expectedWeightsTest4[4] = 0.2158726906203832;
    expectedWeightsTest4[5] = 0.1384130236945297;
    expectedWeightsTest4[6] = 0.02380952381213262;

    scalarList expectedAbscissaeTest4(7);

    expectedAbscissaeTest4[0] = 0.0;
    expectedAbscissaeTest4[1] = 0.08488805184722899;
    expectedAbscissaeTest4[2] = 0.2655756032333322;
    expectedAbscissaeTest4[3] = 0.499999999964492;
    expectedAbscissaeTest4[4] = 0.734424396710906;
    expectedAbscissaeTest4[5] = 0.9151119481304069;
    expectedAbscissaeTest4[6] = 1.0;

    testQuadrature
    (
        mGaussTest4,
        quadraturePropertiesLobatto,
        expectedWeightsTest4,
        expectedAbscissaeTest4,
        "Gauss-Lobatto"
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
