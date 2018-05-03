/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 Alberto Passalacqua
     \\/     M anipulation  |
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
    Application to construct moments give other parameters. Builds initial
    field files for main solver to use as input.

Description
    Preprocessing application to eliminate the need to create fields for all
    moments. Instead moments are consucted using inputs from
    momentGenerationDict. Different methods can be used.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "univariateMomentSet.H"
#include "extendedMomentInversion.H"
#include "mappedPtrList.H"
#include "probes.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addTimeOptions.H"
    #include "setRootCase.H"
    #include "createTime.H"

    // List of times
    instantList Times = runTime.times();

    // Default to latest time
    label startTime = Times.size() - 1;
    if (args.optionFound("time"))
    {
        Foam::scalar timeValue = args.optionRead<scalar>("time");

        startTime = Foam::Time::findClosestTimeIndex(Times, timeValue);
    }
    runTime.setTime(Times[startTime], startTime);

    #include "createMesh.H"

    Info<< "Reading pointDistributionDict" << nl << endl;
    IOdictionary distDict
    (
        IOobject
        (
            "pointDistributionDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<word> phases(distDict.lookup("phases"));

    for (label phasei = 0; phasei < phases.size(); phasei++)
    {
        word phaseName = phases[phasei];
        const dictionary& phaseDict(distDict.subDict(phaseName));
        const dictionary& probesDict(phaseDict.subDict("probes"));

        Info<< "Reading quadratureProperties." << phaseName << endl;
        IOdictionary quadratureDict
        (
            IOobject
            (
                IOobject::groupName
                (
                    "quadratureProperties",
                    phaseName
                ),
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelListList momentOrders(quadratureDict.lookup("moments"));
        labelListList nodeIndexes(quadratureDict.lookup("nodes"));
        label nMoments = momentOrders.size();
        label nSamples = phaseDict.lookupOrDefault("nSamples", 100);
        label nSecondaryNodes =
            quadratureDict.subDict("extendedMomentInversion"). lookupOrDefault("nSecondaryNodes", nMoments + 1);

        // CreateProbes and probesDict
        fileName probeDir;
        fileName probeSubDir = phaseName + "Probes";
        probes mProbes(probeSubDir, runTime, probesDict);

        if (mesh.name() != polyMesh::defaultRegion)
        {
            probeSubDir = probeSubDir/mesh.name();
        }
        probeSubDir = "postProcessing"/probeSubDir/mesh.time().timeName();
        probeDir = runTime.path()/probeSubDir;
        // Remove ".."
        probeDir.clean();
        mkDir(probeDir);

        unsigned int p = IOstream::defaultPrecision() + 7;

        OFstream outputFile(probeSubDir/"quadrature");
        outputFile  << "# Quadrature" << nl
                    << '#' << setw(p - 1)
                    << "abscissae" << ' ' << setw(p - 9)
                    << "n" << endl;


        // Create moment set where each entry is the list of probed moments
        mappedList<scalarList> momentProbes
        (
            nMoments,
            momentOrders,
            scalarList(mProbes.size())
        );

        // Read moments in from fields and copy probed location
        forAll(momentOrders, mi)
        {
            const labelList& momentOrder = momentOrders[mi];
            volScalarField momenti
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        IOobject::groupName
                        (
                            "moment",
                            mappedPtrList<scalar>::listToWord(momentOrder)
                        ),
                        phaseName
                    ),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            momentProbes(momentOrder) = mProbes.sample(momenti);
        }

        scalarListList weights(mProbes.size());
        scalarListList abscissae(mProbes.size());

        forAll(mProbes, probei)
        {

            autoPtr<extendedMomentInversion> EQMOM
            (
                extendedMomentInversion::New
                (
                    quadratureDict.subDict("extendedMomentInversion"),
                    nMoments,
                    nSecondaryNodes
                )
            );

            univariateMomentSet moments(nMoments, "RPlus");
            forAll(moments, mi)
            {
                moments(momentOrders[mi][0]) =
                    momentProbes(momentOrders[mi])[probei];
            }

            //- Guarentee EQMOM is used by setting an odd number of moments
            label nRealizableMoments = moments.nRealizableMoments();
            if (nRealizableMoments % 2 == 0)
            {
                nRealizableMoments--;
            }
            for (label mi = nRealizableMoments; mi < nMoments; mi++)
            {
                moments(momentOrders[mi][0]) = 0.0;
            }

            //- Invert moments
            EQMOM->invert(moments);

            const scalarRectangularMatrix& sAbscissae
            (
                EQMOM->secondaryAbscissae()
            );

            //- COnstruct sample distribution
            scalarField x(nSamples, 0.0);
            scalar xMax = phaseDict.lookupOrDefault("xMax", max(sAbscissae));
            scalar xMin = phaseDict.lookupOrDefault("xMin", min(sAbscissae));
            scalar dx = (xMax - xMin)/(nSamples - 1);

            forAll(x, i)
            {
                x[i] = xMin + dx*i;
            }

            scalarField w(EQMOM->f(x)/moments[labelList(1, 0)[0]]);

            forAll(w, i)
            {
                outputFile  << ' ' << setw(p) << x[i]
                            << ' ' << setw(p) << w[i]
                            << endl;
            }
        }
    }

    Info<< nl << "End\n" << endl;

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
