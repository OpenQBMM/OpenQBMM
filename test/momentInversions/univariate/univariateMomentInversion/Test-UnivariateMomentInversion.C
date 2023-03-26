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
#include "newUnivariateMomentInversion.C"

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
    string quadratureName,
    label nMaxNodes = 0
)
{
    showInputMoments(m, quadratureName);

    autoPtr<univariateMomentInversion> inversion
    (
        univariateMomentInversion::New(dict, nMaxNodes)
    );

    inversion().invert(m, 0, 1);

    scalarList weights(inversion().weights());
    scalarList abscissae(inversion().abscissae());

    Info << "Weights: " << weights << endl;
    Info << "Abscissae: " << abscissae << endl;

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

    dictionary quadraturePropertiesGQMOM
    (
        IFstream("quadraturePropertiesGQMOM")()
    );

    // Test 1 - m = (1, 1, 1, 1)
    scalarList inputMoments1(4, 1.0);
    univariateMomentSet mGaussTest1(inputMoments1, "RPlus", SMALL, SMALL);
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

    univariateMomentSet mGaussTest2(inputMoments2, "RPlus", SMALL, SMALL);
    
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

    univariateMomentSet mGaussTest3(inputMoments3, "RPlus", SMALL, SMALL, 1);
    
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

    univariateMomentSet mGaussTest4(inputMoments4, "RPlus", SMALL, SMALL, 2);
    
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

    // Test 5 - GQMOM on R
    scalarList inputMoments5(10);

    inputMoments5[0] = 1.0;
    inputMoments5[1] = 0.0;
    inputMoments5[2] = 1.0;
    inputMoments5[3] = 0.0;
    inputMoments5[4] = 3.0;
    inputMoments5[5] = 0.0;
    inputMoments5[6] = 15.0;
    inputMoments5[7] = 0.0;
    inputMoments5[8] = 105.0;
    inputMoments5[9] = 0.0;
    //inputMoments5[10] = 945.0;
    
    univariateMomentSet mGaussTest5(inputMoments5, "R", SMALL, SMALL, 5);
    
    scalarList expectedWeightsTest5(10, 0);

    expectedWeightsTest5[0] = 4.310652630718267e-06; 
    expectedWeightsTest5[1] = 0.0007580709343122131; 
    expectedWeightsTest5[2] = 0.01911158050077029; 
    expectedWeightsTest5[3] = 0.1354837029802678;
    expectedWeightsTest5[4] = 0.3446423349320191;
    expectedWeightsTest5[5] = 0.3446423349320191;
    expectedWeightsTest5[6] = 0.1354837029802678;
    expectedWeightsTest5[7] = 0.01911158050077032;
    expectedWeightsTest5[8] = 0.0007580709343122137;
    expectedWeightsTest5[9] = 4.310652630718305e-06;

    scalarList expectedAbscissaeTest5(10, 0);

    expectedAbscissaeTest5[0] = -4.859462828332314;
    expectedAbscissaeTest5[1] = -3.581823483551925;
    expectedAbscissaeTest5[2] = -2.484325841638955;
    expectedAbscissaeTest5[3] = -1.465989094391158;
    expectedAbscissaeTest5[4] = -0.4849357075154974;
    expectedAbscissaeTest5[5] = 0.4849357075154979;
    expectedAbscissaeTest5[6] = 1.465989094391158;
    expectedAbscissaeTest5[7] = 2.484325841638951;
    expectedAbscissaeTest5[8] = 3.581823483551929;
    expectedAbscissaeTest5[9] = 4.85946282833231;
    
    testQuadrature
    (
        mGaussTest5,
        quadraturePropertiesGQMOM,
        expectedWeightsTest5,
        expectedAbscissaeTest5,
        "GQMOM",
        10
    );
    

    // Test 6 - GQMOM on R+
    scalarList inputMoments6(10);

    for (label mi = 1; mi < inputMoments6.size() + 1; mi++)
    {
        inputMoments6[mi - 1] = 1.0/scalar(mi);
    }

    univariateMomentSet mGaussTest6(inputMoments6, "RPlus", SMALL, SMALL, 5);
    
    scalarList expectedWeightsTest6(10, 0);

    expectedWeightsTest6[0] = 0.05425312759509376;
    expectedWeightsTest6[1] = 0.0831424672400455;
    expectedWeightsTest6[2] = 0.1350702286390351;
    expectedWeightsTest6[3] = 0.151991920210894;
    expectedWeightsTest6[4] = 0.1773272901966484;
    expectedWeightsTest6[5] = 0.1619822708924677;
    expectedWeightsTest6[6] = 0.1367233648357529;
    expectedWeightsTest6[7] = 0.09873055145020466;
    expectedWeightsTest6[8] = 0.0007787315954269877;
    expectedWeightsTest6[9] = 4.734443075952662e-08;

    scalarList expectedAbscissaeTest6(10, 0);

    expectedAbscissaeTest6[0] = 0.02276445906247883;
    expectedAbscissaeTest6[1] = 0.09312626860812966;
    expectedAbscissaeTest6[2] = 0.2012821161975949;
    expectedAbscissaeTest6[3] = 0.3476412348445122;
    expectedAbscissaeTest6[4] = 0.511307436092814;
    expectedAbscissaeTest6[5] = 0.6864568463904339;
    expectedAbscissaeTest6[6] = 0.8326551774279829;
    expectedAbscissaeTest6[7] = 0.9566638853707417;
    expectedAbscissaeTest6[8] = 1.076481331174686;
    expectedAbscissaeTest6[9] = 1.422414895629834;

    testQuadrature
    (
        mGaussTest6,
        quadraturePropertiesGQMOM,
        expectedWeightsTest6,
        expectedAbscissaeTest6,
        "GQMOM",
        10
    );

    // Test 7 - GQMOM on [0, 1]
    scalarList inputMoments7(10);

    inputMoments7[0] = 1.0/2.0;
    inputMoments7[1] = 3.0/10.0;
    inputMoments7[2] = 1.0/5.0;
    inputMoments7[3] = 1.0/7.0;
    inputMoments7[4] = 3.0/28.0;
    inputMoments7[5] = 1.0/12.0;
    inputMoments7[6] = 1.0/15.0;
    inputMoments7[7] = 3.0/55.0;
    inputMoments7[8] = 1.0/22.0;
    inputMoments7[9] = 1.0/26.0;

    univariateMomentSet mGaussTest7(inputMoments7, "01", SMALL, SMALL, 5);
    
    scalarList expectedWeightsTest7(10, 0);
    
    expectedWeightsTest7[0] = 0.001943873229307942;
    expectedWeightsTest7[1] = 0.007902072129949747; 
    expectedWeightsTest7[2] = 0.02732546156051659;
    expectedWeightsTest7[3] = 0.05761179944540628;
    expectedWeightsTest7[4] = 0.09028841825500214;
    expectedWeightsTest7[5] = 0.10986807966067;
    expectedWeightsTest7[6] = 0.09845777421751999;
    expectedWeightsTest7[7] = 0.06927178902222311;
    expectedWeightsTest7[8] = 0.02951883479525387;
    expectedWeightsTest7[9] = 0.007811897684150718;

    scalarList expectedAbscissaeTest7(10, 0);
    
    expectedAbscissaeTest7[0] = 0.06964594081013871;
    expectedAbscissaeTest7[1] = 0.14017230848925;
    expectedAbscissaeTest7[2] = 0.2355328726617618;
    expectedAbscissaeTest7[3] = 0.3490332042621505;
    expectedAbscissaeTest7[4] = 0.4750293831698755;
    expectedAbscissaeTest7[5] = 0.6015400657359764;
    expectedAbscissaeTest7[6] = 0.723885477138646;
    expectedAbscissaeTest7[7] = 0.8284968176504827;
    expectedAbscissaeTest7[8] = 0.9147863479347316;
    expectedAbscissaeTest7[9] = 0.9684906429192973;

    testQuadrature
    (
        mGaussTest7,
        quadraturePropertiesGQMOM,
        expectedWeightsTest7,
        expectedAbscissaeTest7,
        "GQMOM",
        10
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
