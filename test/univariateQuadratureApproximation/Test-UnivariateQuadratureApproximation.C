/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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
    Test-UnivariateQuadratureApproximation

Description
    Test univariateQuadratureApproximation class and methods.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    label nPrimaryNodes = quadrature.nodes().size();
    label nSecondaryNodes = quadrature.nodes()[0].nSecondaryNodes();

    Info << "\nNumber of primary nodes: " << nPrimaryNodes << endl;

    Info << "\nNumber of secondary nodes: " << nSecondaryNodes << endl;

    Info << "\nInverting moments. " << endl;

    quadrature.updateQuadrature();

    Info<< "\nStoring quadrature fields.\n" << endl;

    for (label nodeI = 0; nodeI < nPrimaryNodes; nodeI++)
    {
        quadrature.nodes()[nodeI].primaryWeight().write();
        quadrature.nodes()[nodeI].primaryAbscissa().write();
        quadrature.nodes()[nodeI].sigma().write();

        for (label sNodeI = 0; sNodeI < nSecondaryNodes; sNodeI++)
        {
            quadrature.nodes()[nodeI].secondaryWeights()[sNodeI].write();
            quadrature.nodes()[nodeI].secondaryAbscissae()[sNodeI].write();
        }
    }

    runTime++;

    quadrature.updateMoments();

    for (label mI = 0; mI < quadrature.nMoments(); mI++)
    {
        quadrature.moments()[mI].write();
    }

    runTime++;

    for (label mI = 0; mI < quadrature.nMoments(); mI++)
    {
        quadrature.moments()[mI] *= 2.0;
    }

    quadrature.updateQuadrature();
    quadrature.updateMoments();

    for (label mI = 0; mI < quadrature.nMoments(); mI++)
    {
        quadrature.moments()[mI].write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
