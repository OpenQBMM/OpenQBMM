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
    Test-GaussMomentInversion

Description
    Test gaussMomentInversion class and methods.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "scalarMatrices.H"
#include "univariateMomentSet.H"
#include "gaussMomentInversion.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info << "Testing GaussMomentInversion\n" << endl;

    label nMoments = 4;

    labelListList momentOrders(0);

    for(label mi = 0; mi < nMoments; mi++)
    {
        labelList a(1, mi);

        momentOrders.append(a);
    }

    Info << "List of moment orders" << endl << momentOrders << endl;

    univariateMomentSet m
    (
        nMoments,
        momentOrders,
        "R"
    );

    Info << "Moment set (empty): " << m << endl;

    m[0] = 1.0;
    m[1] = 2.708217669;
    m[2] = 8.951330468;
    m[3] = 35.95258119;

    Info << setprecision(16);
    Info << "\nInput moments\n" << endl;

    for (label momenti = 0; momenti < nMoments; momenti++)
    {
        Info << "Moment " << momenti << " = " << m[momenti] << endl;
    }

    gaussMomentInversion quadrature(m);
    quadrature.invert();

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

    Info << "The number of invertible moments is "
         << quadrature.nInvertibleMoments() << "\n" << endl;

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
