/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 Alberto Passalacqua
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    reconstructPointDistribution

Description
    Utility to computes moments.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Utility to computes moments."
    );

    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // List of times
    instantList Times = runTime.times();

    // Default to latest time
    label startTime = Times.size() - 1;

    if (args.found("time"))
    {
        Foam::scalar timeValue = args.get<scalar>("time");
        startTime = Foam::Time::findClosestTimeIndex(Times, timeValue);
    }

    runTime.setTime(Times[startTime], startTime);

    #include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "computeMomentsDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    wordList phaseNames(dict.lookup("phases"));

    forAll(phaseNames, phasei)
    {
        velocityQuadratureApproximation quadrature
        (
            phaseNames[phasei],
            mesh,
            "RPlus"
        );

        autoPtr<mappedPtrList<volVelocityNode>> nodes(&(quadrature.nodes()));

        PtrList<dimensionSet> abscissaeDimensions
        (
            quadrature.momentOrders()[0].size()
        );

        labelList scalarIndexes = nodes()[0].scalarIndexes();

        if (scalarIndexes.size() == 0)
        {
            scalarIndexes.append(-1);
        }

        label si = 0;

        forAll(quadrature.momentOrders()[0], cmpti)
        {
            if (cmpti == scalarIndexes[si])
            {
                abscissaeDimensions.set
                (
                    cmpti,
                    new dimensionSet(nodes()[0].primaryAbscissae()[si].dimensions())
                );

                si++;
            }
            else
            {
                abscissaeDimensions.set
                (
                    cmpti,
                    new dimensionSet(dimVelocity)
                );
            }
        }

        labelListList newMomentOrders
        (
            dict.subDict(phaseNames[phasei]).lookup("moments")
        );

        forAll(newMomentOrders, mi)
        {
            dimensionSet mDims(nodes()[0].primaryWeight().dimensions());

            forAll(abscissaeDimensions, cmpti)
            {
                mDims *= pow(abscissaeDimensions[cmpti], newMomentOrders[mi][cmpti]);
            }

            volVelocityMoment moment
            (
                phaseNames[phasei],
                newMomentOrders[mi],
                nodes,
                volScalarField
                (
                    IOobject
                    (
                        "tmp",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar("0", mDims, 0.0),
                    quadrature.moments()[0].boundaryField().types()
                )
            );

            Info<< "Created " << moment.name() << endl;

            moment.update();
            moment.updateBoundaries();
            moment.write();
        }
    }

    Info<< nl << "End\n" << endl;

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
