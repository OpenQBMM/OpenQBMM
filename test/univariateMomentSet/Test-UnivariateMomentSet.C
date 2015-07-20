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

#include "IOmanip.H"
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
	
// 	m[4] = -10.0;
// 	m[5] = 0.000000001;
// 	m[6] = 0.111111111;

	Info << setprecision(16);
    Info << "Input moments\n" << endl;

	for (label momentI = 0; momentI < nMoments; momentI++)
	{
	    Info << "Moment " << momentI << " = " << m[momentI] << endl; 
	}

	univariateMomentSet moments(m);
    moments.invert();
	
	if(moments.isRealizable())
	{
		Info << "\nThe full set of moments is realizable.\n" << endl ;
	}
	else
	{
		Info << "\nThe full set of moments is not realizable.\n" << endl
             << "The number of realizable moments is "
             << moments.nRealizableMoments() << "\n" << endl;
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
