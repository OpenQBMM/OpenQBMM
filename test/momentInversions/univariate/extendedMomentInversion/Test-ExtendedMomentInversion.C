/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2026 Alberto Passalacqua
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
    Test-ExtendedMomentInversion

Description
    Test the extendedMomentInversion class and its kernel density functions.

    The check that carries the weight is the reconstruction of the moments:
    the extended quadrature method of moments determines n weights, n
    abscissae and sigma from 2 n + 1 moments, so the quadrature it builds has
    to reproduce every one of them. That is a property of the method, not a
    value recorded from a previous run, and it fails for any corruption of
    the primary quadrature, of the parameter of the kernel, or of the
    secondary quadrature.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "IFstream.H"
#include "scalarList.H"
#include "scalarMatrices.H"
#include "supportType.H"
#include "univariateMomentSet.H"
#include "extendedMomentInversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void compareScalar
(
    const scalar computed,
    const scalar expected,
    const scalar tolerance,
    const string& name
)
{
    Info<< "  " << name << " = " << computed
        << ", expected " << expected << endl;

    const scalar difference = mag(computed - expected);

    if (difference > tolerance*max(mag(expected), SMALL))
    {
        FatalErrorInFunction
            << "The value of " << name << " is not the expected one." << nl
            << "    Computed: " << computed << nl
            << "    Expected: " << expected << nl
            << "    Difference: " << difference << nl
            << "    Tolerance: " << tolerance << nl
            << exit(FatalError);
    }
}


//- Moment of order momentOrder of the reconstructed distribution, computed
//  from the primary and secondary quadrature
Foam::scalar momentFromQuadrature
(
    const extendedMomentInversion& EQMOM,
    const label momentOrder
)
{
    const scalarList& pWeights(EQMOM.primaryWeights());
    const scalarRectangularMatrix& sWeights(EQMOM.secondaryWeights());
    const scalarRectangularMatrix& sAbscissae(EQMOM.secondaryAbscissae());

    scalar moment = 0.0;

    for (label pNodei = 0; pNodei < EQMOM.nPrimaryNodes(); pNodei++)
    {
        scalar secondarySum = 0.0;

        for (label sNodei = 0; sNodei < EQMOM.nSecondaryNodes(); sNodei++)
        {
            secondarySum +=
                sWeights[pNodei][sNodei]
               *pow(sAbscissae[pNodei][sNodei], momentOrder);
        }

        moment += pWeights[pNodei]*secondarySum;
    }

    return moment;
}


void showInputMoments(const scalarList& moments, const string& name)
{
    Info<< "\nTesting " << name << "\n" << endl;

    forAll(moments, momenti)
    {
        Info<< "  inputMoments[" << momenti << "] = " << moments[momenti]
            << endl;
    }

    Info<< endl;
}


//- Invert a moment vector and check that the quadrature reproduces it.
//  nPreservedMoments is the number of moments the method conserves, which is
//  all of them when a valid sigma is found.
void testEQMOM
(
    const dictionary& dict,
    const word& kernel,
    const scalarList& inputMoments,
    const supportType& support,
    const scalar expectedSigma,
    const label nPreservedMoments,
    const scalar tolerance
)
{
    dictionary quadratureDict(dict);
    quadratureDict.set("extendedMomentInversion", kernel);

    showInputMoments(inputMoments, kernel + " kernel density function");

    univariateMomentSet moments(inputMoments, support, SMALL, SMALL);

    autoPtr<extendedMomentInversion> EQMOM
    (
        extendedMomentInversion::New
        (
            quadratureDict,
            inputMoments.size(),
            readLabel(quadratureDict.lookup("nSecondaryNodes"))
        )
    );

    EQMOM->invert(moments);

    Info<< "\nVerifying the parameter of the kernel density function\n" << endl;

    compareScalar(EQMOM->sigma(), expectedSigma, tolerance, "sigma");

    Info<< "\nVerifying moment conservation\n" << endl;

    for (label momenti = 0; momenti < nPreservedMoments; momenti++)
    {
        compareScalar
        (
            momentFromQuadrature(EQMOM(), momenti),
            inputMoments[momenti],
            tolerance,
            "moment " + Foam::name(momenti)
        );
    }

    Info<< endl;
}


//- A moment vector of a Dirac delta has no spread for the kernel density
//  function to represent, so sigma has to be null and the reconstruction has
//  to fall back on the primary quadrature alone
void testDiracDelta
(
    const dictionary& dict,
    const word& kernel,
    const scalarList& inputMoments,
    const supportType& support
)
{
    dictionary quadratureDict(dict);
    quadratureDict.set("extendedMomentInversion", kernel);

    showInputMoments
    (
        inputMoments,
        kernel + " kernel density function, Dirac delta"
    );

    univariateMomentSet moments(inputMoments, support, SMALL, SMALL);

    autoPtr<extendedMomentInversion> EQMOM
    (
        extendedMomentInversion::New
        (
            quadratureDict,
            inputMoments.size(),
            readLabel(quadratureDict.lookup("nSecondaryNodes"))
        )
    );

    EQMOM->invert(moments);

    Info<< "\nVerifying the fall back to the primary quadrature...";

    if (EQMOM->sigma() != 0.0)
    {
        FatalErrorInFunction
            << "The kernel density function of a Dirac delta has a non-null "
            << "parameter." << nl
            << "    sigma: " << EQMOM->sigma() << nl
            << exit(FatalError);
    }

    Info<< "OK" << endl;

    Info<< "\nVerifying moment conservation\n" << endl;

    // The zero and first order moments are the ones a single node of the
    // primary quadrature determines
    for (label momenti = 0; momenti < 2; momenti++)
    {
        compareScalar
        (
            momentFromQuadrature(EQMOM(), momenti),
            inputMoments[momenti],
            1.0e-10,
            "moment " + Foam::name(momenti)
        );
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<< "Reading quadratureProperties\n" << endl;

    IFstream quadraturePropertiesFile("quadratureProperties");
    dictionary quadratureProperties(quadraturePropertiesFile);

    Info<< setprecision(16);

    // Moments of a lognormal distribution of parameter sigma = 0.4, the
    // moment vector this test has always used. The kernel density function
    // of the same family has to recover that parameter, to the accuracy of
    // the ten significant digits the moments are given with.
    const scalarList lognormalMoments
    ({
        1.0, 2.708217669, 8.951330468, 35.95258119, 174.4370267
    });

    testEQMOM
    (
        quadratureProperties,
        "lognormal",
        lognormalMoments,
        supportType::RPlus,
        0.4,
        5,
        1.0e-8
    );

    // A gamma kernel density function on the same moment vector has no
    // closed form to recover, so the value is the one the method converges
    // to. It is here to catch a change of the search for sigma.
    testEQMOM
    (
        quadratureProperties,
        "gamma",
        lognormalMoments,
        supportType::RPlus,
        0.4379509305189075,
        5,
        1.0e-8
    );

    // Exact moments of a beta distribution of shape a = 2 and b = 3, which
    // has support over [0, 1]. The kernel density function of the same
    // family has to recover it exactly, with a parameter
    // sigma = 1/(a + b + 1) = 1/6, and a single primary node carrying it.
    const scalarList betaMoments
    ({
        1.0,
        0.4,
        0.2,
        0.11428571428571428,
        0.07142857142857142
    });

    testEQMOM
    (
        quadratureProperties,
        "beta",
        betaMoments,
        supportType::ZeroOne,
        1.0/6.0,
        5,
        1.0e-12
    );

    // Moments of a Dirac delta at exp(1/2), which the commented moment
    // vectors of this test used to cover
    const scalarList diracMoments
    ({
        1.0, 1.6487212707, 2.7182818285, 4.4816890703, 7.3890560989
    });

    testDiracDelta
    (
        quadratureProperties,
        "lognormal",
        diracMoments,
        supportType::RPlus
    );

    testDiracDelta
    (
        quadratureProperties,
        "gamma",
        diracMoments,
        supportType::RPlus
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
