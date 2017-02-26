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
    Test-GaussLobattoMomentInversion

Description
    Test gaussLobattoMomentInversion class and methods.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "scalarMatrices.H"
#include "univariateMomentSet.H"
#include "gaussLobattoMomentInversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info << "Testing GaussLobattoMomentInversion\n" << endl;

    label nMoments = 6;

    univariateMomentSet m
    (
        nMoments,
        "R",
        2
    );

    Info << "Moment set (empty): " << m << endl;

    m[0] = 1.0000;
    m[1] = 0.2500;
    m[2] = 0.1058;
    m[3] = 0.0562;
    m[4] = 0.0340;
    m[5] = 0.0224;

    Info << setprecision(16);
    Info << "\nInput moments\n" << endl;

    for (label momenti = 0; momenti < nMoments; momenti++)
    {
        Info << "Moment " << momenti << " = " << m[momenti] << endl;
    }

    if (m.isFullyRealizable())
    {
        Info << "\nThe full set of moments is realizable.\n" << endl ;
    }
    else if (m.isSubsetRealizable())
    {
        Info << "\nThe full set of moments is not realizable.\n" << endl
             << "The number of realizable moments is "
             << m.nRealizableMoments() << "\n" << endl;
    }
    else
    {
        Info << "\nThe moment set is not realizable.\n" << endl;
    }

    gaussLobattoMomentInversion quadrature(m, 0, 1);

    Info << "The number of invertible moments is "
         << quadrature.nInvertibleMoments() << "\n" << endl;

    quadrature.invert();

    scalarList weights(quadrature.weights());
    scalarList abscissae(quadrature.abscissae());

    Info << "Weights and abscissae:\n" << endl;

    forAll (weights, nodei)
    {
        Info << "Node " << nodei
             << " Weight: " << weights[nodei]
             << " Abscissa: " << abscissae[nodei] << endl;
    }

    m.update(weights, abscissae);

    Info << "\nMoments computed from quadrature\n" << endl;

    for (label momenti = 0; momenti < nMoments; momenti++)
    {
        Info << "Moment " << momenti << " = " << m[momenti] << endl;
    }

    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
