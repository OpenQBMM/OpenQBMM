/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-07-20 Alberto Passalacqua: Implemented 2D moment initial condition to
                                test moment advection schemes.
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

Description
    Set initial condition for moments in a 2D problem to test moment advection
    schemes.

    The computational domain is assumed to extend along the x and y directions,
    with cell centers having 0 < x < 1; 0 < y < 1.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("setMoments2DDict");
    #include "setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictName << "\n" << endl;

    IOdictionary setMoments1DDict(dictIO);

    label nMoments = setMoments1DDict.lookupOrDefault<label>("nMoments", 4);

    word initType
        = setMoments1DDict.lookupOrDefault<word>("initType", "regular");

    Info << "Number of moments to initialize: " << nMoments << endl;

    PtrList<volScalarField> moments(nMoments);

    forAll(moments, mi)
    {
        moments.set
        (
            mi,
            new volScalarField
            (
                IOobject
                (
                    "moment." + Foam::name(mi) + ".populationBalance",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            )
        );
    }

    volVectorField cellCentres(mesh.C());

    Info << "Initialising moments..." << endl;

    if (initType == "regular")
    {
        forAll(cellCentres, celli)
        {
            scalar x = cellCentres[celli].x();
            scalar y = cellCentres[celli].y();
            scalar twoPi = 2.0*Foam::constant::mathematical::pi;

            scalar z = 8.0*Foam::sqrt(sqr(x - 1.0/8.0) + sqr(y - 1.0/8.0));
            scalar oneMinusZ = 1.0 - z;

            forAll(moments, mi)
            {
                moments[mi][celli] = 0.0;
            }

            if (z >= 0.0 && z <= 1.0)
            {
                scalar alpha = 3.5 + 1.5*Foam::sin(twoPi*oneMinusZ);
                scalar beta = 3.5 - 1.5*Foam::cos(twoPi*oneMinusZ);
                moments[0][celli] = 1.0;

                for (label mi = 1; mi < nMoments; mi++)
                {
                    moments[mi][celli] =
                        moments[mi - 1][celli]*(alpha + scalar(mi - 1))
                    /(alpha + beta + scalar(mi - 1));
                }

                scalar coeff = 0.5*
                (
                    1.0 + Foam::tanh(Foam::tan(Foam::constant::mathematical::pi
                    *(oneMinusZ - 0.5)))
                );

                for (label mi = 0; mi < nMoments; mi++)
                {
                    moments[mi][celli] *= coeff;
                }
            }
        }
    }
    else if (initType == "bimodal")
    {
        NotImplemented;
    }
    else
    {
        Info << "Valid initialisation types: regular, bimodal." << endl;
        exit(1);
    }

    Info << "Writing moments..." << endl;

    forAll(moments, mi)
    {
        moments[mi].write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
