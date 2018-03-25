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
    Solver for a system of any number of compressible fluid phases with a
    common pressure, but otherwise separate properties. The type of phase model
    is run time selectable and can optionally represent multiple species and
    in-phase reactions. The phase system is also run time selectable and can
    optionally represent different types of momentun, heat and mass transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "momentGenerationModel.H"

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

        // Read number of nodes from quadratureProperties.phase
        label nNodes;

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

        nNodes = HashTable<dictionary>(quadratureDict.lookup("nodes")).size();

        label nMoments = 2*nNodes;
        bool radau = phaseDict.lookupOrDefault<bool>("Radau", false);
        bool extended(phaseDict.lookupOrDefault<bool>("extended", false));

        if (radau)
        {
            nNodes++;
        }

        if (extended || radau)
        {
            nMoments++;
        }

        autoPtr<momentGenerationModel> momentGenerator
            = momentGenerationModel::New(phaseDict, nNodes, extended, radau);

        PtrList<volScalarField> moments(nMoments);

        //  Set internal field values and initialize moments.
        {
            const dictionary& dict(phaseDict.subDict("internal"));

            momentGenerator().updateQuadrature(dict);

            forAll(moments, mi)
            {
                moments.set
                (
                    mi,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "moment",
                                IOobject::groupName
                                (
                                    Foam::name(mi),
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

                moments[mi].boundaryFieldRef().readField
                (
                    moments[mi].internalField(),
                    boundaries
                );

                Info<< "    " << moments[mi].name() << endl;
            }
        }

        forAll(mesh.boundaryMesh(), bi)
        {
            word bName = mesh.boundaryMesh()[bi].name();

            if (!phaseDict.found(bName))
            {
                bName = "default";
            }

            if (moments[0].boundaryField()[bi].fixesValue())
            {
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

        forAll(moments, mi)
        {
            moments[mi].write();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}
