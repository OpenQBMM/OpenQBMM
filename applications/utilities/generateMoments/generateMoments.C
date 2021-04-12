/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 Alberto Passalacqua
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    Utility to generate moments to initialize solvers. Instead moments are 
    consucted using inputs from momentGenerationDict. 
    Different methods can be used to specify the moment definition.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "momentGenerationModel.H"
#include "mappedPtrList.H"
#include "topoSetSource.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Utility to generate moments to initialize QBMM solvers."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading initial conditions and creating moments" << nl << endl;

    IOdictionary phaseDicts
    (
        IOobject
        (
            "momentGenerationDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    List<word> phases(phaseDicts.lookup("phases"));

    for (label phasei = 0; phasei < phases.size(); phasei++)
    {
        word phaseName = phases[phasei];
        const dictionary& phaseDict(phaseDicts.subDict(phaseName));

        Info<< "Creating moments for phase: " << phaseName << endl;

        // Read number of nodes from quadratureProperties.phaseName
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
        label nNodes = nodeIndexes.size();

        autoPtr<momentGenerationModel> momentGenerator =
            momentGenerationModel::New
            (
                mesh,
                phaseDict,
                momentOrders,
                nNodes
            );

        mappedPtrList<dictionary> boundaries(nMoments, momentOrders);

        forAll(momentOrders, mi)
        {
            word bName =
                "moment"
              + mappedScalarList::listToWord(momentOrders[mi])
              + "Boundary";

            boundaries.set
            (
                momentOrders[mi],
                new dictionary
                (
                    phaseDict.found(bName)
                  ? phaseDict.subDict(bName)
                  : phaseDicts.subDict("boundaries")
                )
            );
        }


        mappedPtrList<volScalarField> moments(nMoments, momentOrders);

        //  Set internal field values and initialize moments.
        {
            Info<< "Setting internal fields" << nl << endl;

            const dictionary& dict
            (
                phaseDict.found("internal")
              ? phaseDict.subDict("internal")
              : phaseDict.subDict("default")
            );

            momentGenerator().updateMoments(dict);

            forAll(moments, mi)
            {
                const labelList& momentOrder = momentOrders[mi];

                moments.set
                (
                    momentOrder,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "moment",
                                IOobject::groupName
                                (
                                    mappedPtrList<scalar>::listToWord(momentOrder),
                                    phases[phasei]
                                )
                            ),
                            "0",
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "moment",
                            momentGenerator().momentDims()(momentOrder),
                            0.0
                        )
                    )
                );

                moments[mi].primitiveFieldRef() =
                    momentGenerator().moments()(momentOrder);

                //  Set boundaries based oboundary section
                //  Initial values specified in the dictionary are overwritten
                moments(momentOrder).boundaryFieldRef().readField
                (
                    moments(momentOrder).internalField(),
                    boundaries[mi]
                );
            }
        }

        forAll(mesh.boundaryMesh(), bi)
        {
            if
            (
                moments[0].boundaryField()[bi].fixesValue()
             && (
                    phaseDict.found(mesh.boundaryMesh()[bi].name())
                 || phaseDict.found("default")
                )
            )
            {
                Info<< "Setting " << mesh.boundaryMesh()[bi].name()
                    << " boundary" << endl;

                const dictionary& dict
                (
                    phaseDict.found(mesh.boundaryMesh()[bi].name())
                  ? phaseDict.subDict(mesh.boundaryMesh()[bi].name())
                  : phaseDict.subDict("default")
                );

                momentGenerator().updateMoments(dict, bi);

                forAll(moments, mi)
                {
                    moments[mi].boundaryFieldRef()[bi] ==
                        momentGenerator().moments()[mi];
                }
            }
        }

        //- Set regions of domain using methods seen in setFields
        if (phaseDict.found("regions"))
        {
            PtrList<entry> regions(phaseDict.lookup("regions"));

            forAll(regions, regionI)
            {
                const entry& region = regions[regionI];

                autoPtr<topoSetSource> source =
                    topoSetSource::New(region.keyword(), mesh, region.dict());

                if (source().setType() == topoSetSource::CELLSETSOURCE)
                {
                    cellSet selectedCellSet
                    (
                        mesh,
                        "cellSet",
                        mesh.nCells()/10+1  // Reasonable size estimate.
                    );

                    source->applyToSet
                    (
                        topoSetSource::NEW,
                        selectedCellSet
                    );

                    const labelList& cells = selectedCellSet.toc();
                    momentGenerator().updateMoments(region.dict(), cells);

                    forAll(cells, celli)
                    {
                        forAll(moments, mi)
                        {
                            moments[mi][cells[celli]] =
                                momentGenerator().moments()[mi][celli];
                        }
                    }

                }
                else if (source().setType() == topoSetSource::FACESETSOURCE)
                {
                    FatalErrorInFunction
                        << "Moments must be volume fields."
                        << abort(FatalError);
                }
            }
        }

        forAll(moments, mi)
        {
            moments[mi].write();
        }
    }

    Info<< nl << "End\n" << endl;

    return 0;
}
