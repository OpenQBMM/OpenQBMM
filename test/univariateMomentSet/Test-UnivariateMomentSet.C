/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 Alberto Passalacqua
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
    Test-UnivariateMomentSet

Description
    Test univariateMomentSet class and methods.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "scalarMatrices.H"
#include "univariateMomentSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info << "Testing univariateMomentSet\n" << endl;

    label nMoments = 7;

    scalarDiagonalMatrix m(nMoments, 0.0);
    // Computing integer moments of a log-normal function
//     scalarDiagonalMatrix m(nMoments, 1.0);
//     scalar mu = 0.0;
//     scalar sigma = 0.25;
//
//     for (label momentI = 0; momentI < nMoments; momentI++)
//     {
//         m[momentI] = Foam::exp(momentI*mu + Foam::sqr(momentI*sigma)/2.0);
//     }

    // Computing integer moments of a gaussian function
//     scalar mu = 0.0;
//     scalar sigma = 0.25;


//     m[0] = 1.0;
//
//     for (label momentI = 2; momentI < nMoments; momentI = momentI + 2)
//     {
//         m[momentI] = pow(sigma, momentI)*Foam::factorial(momentI-1);
//     }

    // Computing integer moments of sum of two beta function
//     scalar mu1 = 0.5;
//     scalar sigma1 = 0.3;
//     scalar mu2 = 0.5;
//     scalar sigma2 = 0.3;
//
//     m[0] = 1.0;
//
//     for (label momentI = 1; momentI < nMoments; momentI++)
//     {
//         m[momentI] = (((mu1 + (momentI - 1.0)*sigma1)*m[momentI-1]
//             /(1.0 + (momentI - 1.0)*sigma1))
//             +((mu2 + (momentI - 1.0)*sigma2)*m[momentI-1]
//             /(1.0 + (momentI - 1.0)*sigma2)))/2;
//     }

// Computing integer moments of a beta function
//     scalar mu = 0.5;
//     scalar sigma = 0.3;
//
//     m[0] = 1.0;
//
//     for (label momentI = 1; momentI < nMoments; momentI++)
//     {
//         m[momentI] = (mu + (momentI - 1.0)*sigma)*m[momentI-1]
//             /(1.0 + (momentI - 1.0)*sigma);
//     }

//    m[0] = 1.0;
//    m[1] = 0.0;

    m[0] = 3.125e12;
    m[1] = 6.25e6;
    m[2] = 12.5;
    m[3] = 2.5e-5;
    m[4] = 5.0e-11;
    m[5] = 1.0e-16;
    m[6] = 2.0e-22;

    word support = "RPlus";

    Info << "Support: " << support << endl;

    Info << setprecision(16);
    Info << "\nInput moments\n" << endl;

    for (label momentI = 0; momentI < nMoments; momentI++)
    {
        Info << "Moment " << momentI << " = " << m[momentI] << endl;
    }

    univariateMomentSet moments(m, support);

    Info << "\nStored moments\n" << endl;

    forAll(moments, momentI)
    {
        Info << "Moment " << momentI << " = " << moments[momentI] << endl;
    }

    moments.invert();

    if (moments.isFullyRealizable())
    {
	Info << "\nThe full set of moments is realizable.\n" << endl ;
    }
    else if (moments.isSubsetRealizable())
    {
        Info << "\nThe full set of moments is not realizable.\n" << endl
             << "The number of realizable moments is "
             << moments.nRealizableMoments() << "\n" << endl;
    }
    else
    {
        Info << "\nThe moment set is not realizable.\n" << endl;
    }

    Info << "The number of invertible moments is "
	 << moments.nInvertibleMoments() << "\n" << endl;

    scalarDiagonalMatrix weights(moments.weights());
    scalarDiagonalMatrix abscissae(moments.abscissae());

    Info << "Weights and abscissae:\n" << endl;

    for (label nodeI = 0; nodeI < moments.nNodes(); nodeI++)
    {
	Info << "Node " << nodeI
             << " Weight: " << weights[nodeI]
             << " Abscissa: " << abscissae[nodeI] << endl;
    }

    moments.update();

    Info << "\nMoments computed from quadrature\n" << endl;

    for (label momentI = 0; momentI < nMoments; momentI++)
    {
        Info << "Moment " << momentI << " = " << moments[momentI] << endl;
    }

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
