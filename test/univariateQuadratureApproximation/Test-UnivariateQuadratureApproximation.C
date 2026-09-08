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
    Copyright (C) 2019-2025 Alberto Passalacqua
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
#include "supportType.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nTested = 0;

void check(const bool condition, const string& what)
{
    nTested++;

    if (!condition)
    {
        FatalErrorInFunction
            << "Failed: " << what << nl
            << exit(FatalError);
    }

    Info<< "  OK: " << what << endl;
}


//- Compare two fields cell by cell, to a tolerance relative to the magnitude
//  of the expected value
void checkField
(
    const volScalarField& computed,
    const scalarField& expected,
    const scalar tolerance,
    const string& what
)
{
    nTested++;

    scalar worst = 0;
    label worstCell = -1;

    forAll(expected, celli)
    {
        const scalar error =
            mag(computed[celli] - expected[celli])
           /max(mag(expected[celli]), SMALL);

        if (error > worst)
        {
            worst = error;
            worstCell = celli;
        }
    }

    if (worst > tolerance)
    {
        FatalErrorInFunction
            << "Failed: " << what << nl
            << "    Largest relative error " << worst
            << " in cell " << worstCell << nl
            << "    Computed " << computed[worstCell] << nl
            << "    Expected " << expected[worstCell] << nl
            << exit(FatalError);
    }

    Info<< "  OK: " << what
        << " (largest relative error " << worst << ")" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    label nPrimaryNodes = quadrature.nodes().size();

    Info << "\nNumber of primary nodes: " << nPrimaryNodes << endl;
    Info << "\nInverting moments. " << endl;

    // Keep the moments the case was given, so that what the inversion and
    // the reconstruction do to them can be checked rather than only printed
    PtrList<scalarField> initialMoments(quadrature.nMoments());

    for (label mI = 0; mI < quadrature.nMoments(); mI++)
    {
        initialMoments.set
        (
            mI,
            new scalarField(quadrature.moments()[mI].primitiveField())
        );
    }

    quadrature.updateQuadrature();

    Info<< "\nStoring quadrature fields.\n" << endl;

    for (label nodeI = 0; nodeI < nPrimaryNodes; nodeI++)
    {
        quadrature.nodes()[nodeI].weight().write();
        quadrature.nodes()[nodeI].abscissae()[0].write();
    }

    // The quadrature has as many nodes as the moment set can support, so
    // going to the quadrature and back has to return the moments unchanged
    Info<< "\nMoments recovered from the quadrature" << endl;

    runTime++;

    quadrature.updateMoments();

    for (label mI = 0; mI < quadrature.nMoments(); mI++)
    {
        quadrature.moments()[mI].write();

        checkField
        (
            quadrature.moments()[mI],
            initialMoments[mI],
            1.0e-10,
            "moment of order " + Foam::name(mI) + " is recovered"
        );
    }

    // Doubling every moment doubles the measure without moving it, so the
    // abscissae have to come back unchanged and the weights doubled
    Info<< "\nDoubling the measure" << endl;

    PtrList<scalarField> singleWeights(nPrimaryNodes);
    PtrList<scalarField> singleAbscissae(nPrimaryNodes);

    for (label nodeI = 0; nodeI < nPrimaryNodes; nodeI++)
    {
        singleWeights.set
        (
            nodeI,
            new scalarField
            (
                quadrature.nodes()[nodeI].weight().primitiveField()
            )
        );

        singleAbscissae.set
        (
            nodeI,
            new scalarField
            (
                quadrature.nodes()[nodeI].abscissae()[0].primitiveField()
            )
        );
    }

    runTime++;

    for (label mI = 0; mI < quadrature.nMoments(); mI++)
    {
        quadrature.moments()[mI] *= 2.0;
    }

    quadrature.updateQuadrature();
    quadrature.updateMoments();

    for (label nodeI = 0; nodeI < nPrimaryNodes; nodeI++)
    {
        checkField
        (
            quadrature.nodes()[nodeI].abscissae()[0],
            singleAbscissae[nodeI],
            1.0e-10,
            "abscissa of node " + Foam::name(nodeI) + " is unchanged"
        );

        checkField
        (
            quadrature.nodes()[nodeI].weight(),
            2.0*singleWeights[nodeI],
            1.0e-10,
            "weight of node " + Foam::name(nodeI) + " is doubled"
        );
    }

    for (label mI = 0; mI < quadrature.nMoments(); mI++)
    {
        quadrature.moments()[mI].write();

        checkField
        (
            quadrature.moments()[mI],
            2.0*initialMoments[mI],
            1.0e-10,
            "doubled moment of order " + Foam::name(mI) + " is recovered"
        );
    }

    Info<< "\n" << nTested << " checks passed.\n" << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
