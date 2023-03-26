/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 Alberto Passalacqua
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
    Test-UnivariateMomentSet

Description
    Test univariateMomentSet class and methods.

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

void compareZetas
(
    const scalarList& computedZetas, 
    const scalarList& expectedZetas
)
{
    Info << "\nExpected zetas\n" << endl;
    forAll (expectedZetas, zetai)
    {
        Info<< "  expectedZetas[" << zetai << "] = " << expectedZetas[zetai] 
            << ", computedZetas[" << zetai << "] = " << computedZetas[zetai] 
            << endl;
    }

    if (computedZetas.size() != expectedZetas.size())
    {
        FatalErrorInFunction
            << "Zeta vectors have different size: " 
            << endl
            << "  Size of computed zetas vector: " << computedZetas.size()
            << endl
            << "  Size of expected zetas vector: " << expectedZetas.size()
            << endl
            << exit(FatalError);
    }

    forAll (computedZetas, zetai)
    {
        scalar magDiff = mag(computedZetas[zetai] - expectedZetas[zetai]);

        if (magDiff >= SMALL)
        {
            FatalErrorInFunction
                << "Values of zetas do not match: " 
                << endl
                << "  Position: " << zetai
                << endl
                << "  Expected zeta value " << expectedZetas[zetai]
                << endl
                << "  Computed zeta value " << computedZetas[zetai]
                << endl
                << "  Magnitude of the difference " << magDiff
                << endl
                << exit(FatalError);
        }
    }

    Info << "\nZeta values match.\n" << endl;
}

// Testing with the moment vector
//
// m = (1 1 1 1 1 1 1 1 1 1)
//
// The moment vector is degenerate. The realizability test should 
void testUnitMomentVectorRPlus()
{
    const int nMoments = 10;

    // All input moments are 1
    scalarList inputMoments(nMoments);

    Info << "Testing with unit moment vector\n" << endl
         << "-------------------------------" << endl;

    forAll (inputMoments, mi)
    {
        inputMoments[mi] = 1.0;
        Info << "  inputMoments[" << mi << "] = " << inputMoments[mi] << endl;
    }

    // Expected values of zeta_k
    const int nZetas = nMoments - 1;
    scalarList expectedZetas(nZetas, scalar(0));
    expectedZetas[0] = 1.0;

    Info << "\nAllocating univariateMomentSet\n" << endl;
    univariateMomentSet moments(inputMoments, "RPlus", SMALL, SMALL);

    forAll(moments, mi)
    {
        Info << "  moments[" << mi << "] = " << moments[mi] << endl;
    }

    Info << "\nAllocated univariateMomentSet\n" << endl;

    Info << "\nVerifying full realizability...";
    if (moments.isFullyRealizable())
    {
        FatalErrorInFunction
            << "The input moment vector is incorrectly identified as fully "
            << endl
            << "realizable."
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    Info << "Verifying realizability of subset...";
    if (!moments.isSubsetRealizable())
    {
        FatalErrorInFunction
            << "A subset of the input moment vector is not identified as " 
            << "realizable. " 
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    // Test values of zetas
    scalarList computedZetas(moments.zetas());
    compareZetas(computedZetas, expectedZetas);
}

// This function tests a fully realizable moment vector built as
// 
// m_{i-1} = 1/m_i, with i = 1, ..., nMoments - 1
//
void testFullyRealizableMomentVectorRPlus()
{
    const int nMoments = 20;

    // All input moments are 1
    scalarList inputMoments(nMoments);

    Info << "\nTesting with fully realizable moment vector\n" << endl
         << "-------------------------------------------" << endl;

    for (label mi = 1; mi < nMoments + 1; mi++)
    {
        inputMoments[mi - 1] = 1.0/scalar(mi);
        Info << "  inputMoments[" << mi - 1 << "] = " 
             << inputMoments[mi - 1] << endl;
    }

    // Expected values of zeta_k
    const int nZetas = nMoments - 1;
    scalarList expectedZetas(nZetas);
    
    expectedZetas[0] = 0.5;
    expectedZetas[1] = 0.1666666666666666;
    expectedZetas[2] = 0.3333333333333338;
    expectedZetas[3] = 0.199999999999999;
    expectedZetas[4] = 0.2999999999999967;
    expectedZetas[5] = 0.2142857142857517;
    expectedZetas[6] = 0.2857142857140651;
    expectedZetas[7] = 0.2222222222227857;
    expectedZetas[8] = 0.2777777777779511;
    expectedZetas[9] = 0.2272727272565156;
    expectedZetas[10] = 0.2727272728794396;
    expectedZetas[11] = 0.2307692299916414;
    expectedZetas[12] = 0.2692307724113639;
    expectedZetas[13] = 0.2333333295053653;
    expectedZetas[14] = 0.2666666090668923;
    expectedZetas[15] = 0.2352947518514063;
    expectedZetas[16] = 0.2647010523716313;
    expectedZetas[17] = 0.2368683776418143;
    expectedZetas[18] = 0.2630193512153205;

    Info << "\nAllocating univariateMomentSet\n" << endl;
    univariateMomentSet moments(inputMoments, "RPlus", SMALL, SMALL);

    forAll(moments, mi)
    {
        Info << "  moments[" << mi << "] = " << moments[mi] << endl;
    }

    Info << "\nAllocated univariateMomentSet\n" << endl;

    Info << "\nVerifying full realizability...";
    if (!moments.isFullyRealizable())
    {
        FatalErrorInFunction
            << "The moment vector is fully realizable but not identified as "
            << endl
            << "such."
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    Info << "Verifying realizability of subset...";
    if (!moments.isSubsetRealizable())
    {
        FatalErrorInFunction
            << "A subset of the input moment vector is not identified as " 
            << "realizable. " 
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    // Test values of zetas
    scalarList computedZetas(moments.zetas());
    compareZetas(computedZetas, expectedZetas);
}

// This function tests a fully realizable moment vector built as
// 
// m_{i-1} = 1/m_i, with i = 1, ..., nMoments - 1
// m_1 = 0
//
void testSubsetRealizableMomentVectorRPlus()
{
    const int nMoments = 10;

    // All input moments are 1
    scalarList inputMoments(nMoments);

    Info << "\nTesting with moment vector with three realizable moments\n" 
         << endl
         << "---------------------------------------------------------" 
         << endl;

    for (label mi = 1; mi < nMoments + 1; mi++)
    {
        inputMoments[mi - 1] = 1.0/scalar(mi);

        if (mi - 1 == 3)
        {
            inputMoments[mi - 1] = 0.0;
        }

        Info << "  inputMoments[" << mi - 1 << "] = " 
             << inputMoments[mi - 1] << endl;
    }

    // Expected values of zeta_k
    const int nZetas = nMoments - 1;
    scalarList expectedZetas(nZetas, 0.0);
    
    expectedZetas[0] = 0.5;
    expectedZetas[1] = 0.1666666666666666;
    expectedZetas[2] = -2.666666666666667;

    Info << "\nAllocating univariateMomentSet\n" << endl;
    univariateMomentSet moments(inputMoments, "RPlus", SMALL, SMALL);

    forAll(moments, mi)
    {
        Info << "  moments[" << mi << "] = " << moments[mi] << endl;
    }

    Info << "\nAllocated univariateMomentSet\n" << endl;

    Info << "\nVerifying full realizability...";
    if (moments.isFullyRealizable())
    {
        FatalErrorInFunction
            << "The input moment vector is incorrectly identified as fully "
            << endl
            << "realizable."
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    Info << "Verifying realizability of subset...";
    if (!moments.isSubsetRealizable())
    {
        FatalErrorInFunction
            << "A subset of the input moment vector is not identified as " 
            << "realizable. " 
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    // Test values of zetas
    scalarList computedZetas(moments.zetas());
    compareZetas(computedZetas, expectedZetas);
}

// Testing with a fully realizable set of Gaussian moments
void testGaussianMomentsR()
{
    const int nMoments = 9;

    const scalar mu = -2;
    const scalar sigma = 1;

    // All input moments are 1
    scalarList inputMoments(nMoments);

    Info << "\nTesting with Gaussian moments\n" 
         << endl
         << "------------------------------" 
         << endl;

    inputMoments[0] = 1.0;
    inputMoments[1] = mu;
    inputMoments[2] = sqr(mu) + sqr(sigma);
    inputMoments[3] = pow3(mu) + 3.0*mu*sqr(sigma);
    inputMoments[4] = pow4(mu) + 6.0*sqr(mu)*sqr(sigma) + 3.0*pow3(sigma);
    inputMoments[5] = pow5(mu) + 10.0*pow3(mu)*sqr(sigma) + 15.0*mu*pow4(sigma);

    inputMoments[6] = pow6(mu) + 15.0*pow4(mu)*sqr(sigma) 
        + 45.0*sqr(mu)*pow4(sigma) + 15.0*pow6(sigma);

    inputMoments[7] = pow(mu, 7) + 21.0*pow5(mu)*sqr(sigma) 
        + 105.0*pow3(mu)*pow4(sigma) + 105.0*mu*pow6(sigma);

    inputMoments[8] = pow(mu, 8) + 28.0*pow6(mu)*sqr(sigma) 
        + 210.0*pow4(mu)*pow4(sigma) + 420.0*sqr(mu)*pow6(sigma) 
        + 105.0*pow(sigma, 8);
    
    for (label mi = 1; mi < nMoments + 1; mi++)
    {
        Info << "  inputMoments[" << mi - 1 << "] = " 
             << inputMoments[mi - 1] << endl;
    }

    Info << "\nAllocating univariateMomentSet\n" << endl;
    univariateMomentSet moments(inputMoments, "R", SMALL, SMALL);

    forAll(moments, mi)
    {
        Info << "  moments[" << mi << "] = " << moments[mi] << endl;
    }

    Info << "\nAllocated univariateMomentSet\n" << endl;

    Info << "\nVerifying full realizability...";
    if (!moments.isFullyRealizable())
    {
        FatalErrorInFunction
            << "The moment vector is fully realizable but not identified as "
            << endl
            << "such."
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    Info << "Verifying realizability of subset...";
    if (!moments.isSubsetRealizable())
    {
        FatalErrorInFunction
            << "A subset of the input moment vector is not identified as " 
            << "realizable. " 
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }

    Info << "OK" << endl;

    Info << "Zetas not calculated for moments with support on R." << endl;
}

// Testing with a fully realizable set on [0, 1]
void testRealizableCanonicalMoments()
{
    const int nMoments = 9;

    // All input moments are 1
    scalarList inputMoments(nMoments);

    Info << "\nTesting canonical moments\n" 
         << endl
         << "---------------------------" 
         << endl;

    for (label mi = 1; mi < nMoments + 1; mi++)
    {
        inputMoments[mi - 1] = 1.0/scalar(mi);

        Info << "  inputMoments[" << mi - 1 << "] = " 
             << inputMoments[mi - 1] << endl;
    }
    
    Info << "\nAllocating univariateMomentSet\n" << endl;
    univariateMomentSet moments(inputMoments, "01", SMALL, SMALL);

    forAll(moments, mi)
    {
        Info << "  moments[" << mi << "] = " << moments[mi] << endl;
    }

    Info << "\nAllocated univariateMomentSet\n" << endl;

    Info << "\nVerifying full realizability...";
    if (!moments.isFullyRealizable())
    {
        FatalErrorInFunction
            << "The moment vector is fully realizable but not identified as "
            << endl
            << "such."
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    Info << "Verifying realizability of subset...";
    if (!moments.isSubsetRealizable())
    {
        FatalErrorInFunction
            << "A subset of the input moment vector is not identified as " 
            << "realizable. " 
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }

    Info << "OK" << endl;
}

// Testing with a fully realizable set on [0, 1]
void testUnrealizableCanonicalMoments()
{
    const int nMoments = 9;

    // All input moments are 1
    scalarList inputMoments(nMoments);

    Info << "\nTesting canonical moments\n" 
         << endl
         << "---------------------------" 
         << endl;

    for (label mi = 1; mi < nMoments + 1; mi++)
    {
        inputMoments[mi - 1] = 1.0/scalar(mi);

        if (mi - 1 == 2)
        {
            inputMoments[mi - 1] = -1.0;
        }

        Info << "  inputMoments[" << mi - 1 << "] = " 
             << inputMoments[mi - 1] << endl;
    }
    
    Info << "\nAllocating univariateMomentSet\n" << endl;
    univariateMomentSet moments(inputMoments, "01", SMALL, SMALL);

    forAll(moments, mi)
    {
        Info << "  moments[" << mi << "] = " << moments[mi] << endl;
    }

    Info << "\nAllocated univariateMomentSet\n" << endl;

    Info << "\nVerifying full realizability...";
    if (moments.isFullyRealizable())
    {
        FatalErrorInFunction
            << "The input moment vector is incorrectly identified as fully "
            << endl
            << "realizable."
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }
    Info << "OK" << endl;

    Info << "Verifying realizability of subset...";
    if (!moments.isSubsetRealizable())
    {
        FatalErrorInFunction
            << "A subset of the input moment vector is not identified as " 
            << "realizable. " 
            << endl
            << "  Moment vector: " << moments
            << endl
            << "  Zetas: " << moments.zetas()
            << endl
            << exit(FatalError);
    }

    Info << "OK" << endl;
}

int main(int argc, char *argv[])
{
    Info << "Testing univariateMomentSet\n" << endl;
    Info << setprecision(16);

    testGaussianMomentsR();
    testFullyRealizableMomentVectorRPlus();
    testSubsetRealizableMomentVectorRPlus();
    testUnitMomentVectorRPlus();
    testRealizableCanonicalMoments();
    testUnrealizableCanonicalMoments();

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
