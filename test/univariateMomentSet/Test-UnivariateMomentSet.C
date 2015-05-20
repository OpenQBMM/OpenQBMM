/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 Alberto Passalacqua
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

#include "scalarMatrices.H"
#include "univariateMomentSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	
	Info << "Testing momentSet\n" << endl;
	
	label nMoments = 7;
	scalarDiagonalMatrix m(nMoments, 0.0);
	
	// Computing integer moments of a log-normal function
	scalar mu = 0.0;
	scalar sigma = 0.25;
	
	for (label momentI = 0; momentI < nMoments; momentI++)
	{
		m[momentI] = Foam::exp(momentI*mu + Foam::sqr(momentI*sigma)/2.0);
	}
	
	m[0] = 0.567128698550116;
	m[1] = 0.655890848206384;
	m[2] = 0.777476336281080;
	m[3] = 0.947948891168452;
	m[4] = 1.191578422484480;
	
	for (label momentI = 0; momentI < nMoments; momentI++)
	{
	    Info << "Moment " << momentI << " = " << m[momentI] << endl; 
	}

	univariateMomentSet moments(m);
	
	if(moments.isRealizable())
	{
		Info << "\nMoments are realizable.\n" << endl ;
	}
	else
	{
		Info << "Moments are not realizable.\n" << endl;
	}

	label nInvertibleMoments = moments.nInvertibleMoments();
	
	Info << "The number of invertible moments is " 
		 << nInvertibleMoments << "\n" << endl;
	
	// Note that the number of nodes calculated here is the maximum number
	// of nodes allowed. It is not the actual number of nodes that can be
	// found from the nInvertibleMoments. This is done to keep the test
	// general, since this condition may be encountered in practical cases,
	// when a large moment set is not realizable, but one of its subsets is.
	label nNodes = label(nMoments/2);
	
	scalarDiagonalMatrix w(nNodes);
	scalarDiagonalMatrix abs(nNodes);

	moments.invert(w, abs);

	Info << "Weights and abscissae:\n" << endl;
	
	for (label nodeI = 0; nodeI < nNodes; nodeI++)
	{
		Info << "Node " << nodeI 
			 << " Weight: " << w[nodeI] 
			 << " Abscissa: " << abs[nodeI] << endl;
	}
	
    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
