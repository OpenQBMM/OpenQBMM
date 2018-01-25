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
#include "momentGenerationModel.H"
#include "mappedPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
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
    const dictionary& boundaries(phaseDicts.subDict("boundaries"));

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

        autoPtr<momentGenerationModel> momentGenerator
            = momentGenerationModel::New(phaseDict, momentOrders, nNodes);


        mappedPtrList<volScalarField> moments(nMoments, momentOrders);

        //  Set internal field values and initialize moments.
        {
            const dictionary& dict(phaseDict.subDict("internal"));

            momentGenerator().updateQuadrature(dict);

            forAll(moments, mi)
            {
                const labelList& momentOrder= momentOrders[mi];
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
                        momentGenerator().moments()[mi]
                    )
                );

                //  Set boundaries based oboundary section
                //  Initial values specified in the dictionary are overwritten
                moments[mi].boundaryFieldRef().readField
                (
                    moments[mi].internalField(),
                    boundaries
                );
            }
        }

        forAll(mesh.boundaryMesh(), bi)
        {
            if (moments[0].boundaryField()[bi].fixesValue())
            {
                word bName
                (
                    phaseDict.found(mesh.boundaryMesh()[bi].name())
                  ? mesh.boundaryMesh()[bi].name()
                  : "default"
                );
                dictionary dict = phaseDict.subDict(bName);

                momentGenerator().updateQuadrature(dict);

                forAll(moments, mi)
                {
                    forAll(moments[mi].boundaryField()[bi], facei)
                    {
                        moments[mi].boundaryFieldRef()[bi][facei]
                            = (momentGenerator().moments()[mi]).value();
                    }
                }
            }
        }

        Info<< nl << "Writing moments:" << endl;
        forAll(moments, mi)
        {
            Info<< moments[mi].name() << endl;
            moments[mi].write();
        }
    }

    Info<< nl << "End\n" << endl;

    return 0;
}
