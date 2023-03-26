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
        quadrature.nodes()[nodeI].primaryAbscissae()[0].write();
        quadrature.nodes()[nodeI].sigmas()[0].write();

        for (label sNodeI = 0; sNodeI < nSecondaryNodes; sNodeI++)
        {
            quadrature.nodes()[nodeI].secondaryWeights()[0][sNodeI].write();
            quadrature.nodes()[nodeI].secondaryAbscissae()[0][sNodeI].write();
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
