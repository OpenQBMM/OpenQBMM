/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
-------------------------------------------------------------------------------
2015-03-09 Alberto Passalacqua: Implemented 1D moment initial condition to
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
    Set initial condition for moments in a 1D problem to test moment advection
    schemes.

    The computational domain is assumed to extend along the x direction, with
    cell centers having 0 < x < 1.

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

scalar weights2(const scalar x)
{
    if ((3.0*x >= 1.0) && (x <= 1.0))
    {
        return 9.0*sqr(3.0*x - 1.0)*sqr(1.0 - x);
    }

    return scalar(0);
}

scalar parLambda2(const scalar x)
{
    scalar lambda2min = 2.0e-2;
    scalar lambda2max = 0.7;

    if (3.0*x <= 1.0)
    {
        return lambda2min;
    }
    else if (3.0*x <= 2.0)
    {
        return lambda2min*sqr(2.0 - 3.0*x)*(6.0*x - 1.0)
            + lambda2max*sqr(3.0*x - 1.0)*(5.0 - 6.0*x);
    }

    return lambda2max;
}

scalar park2(scalar x)
{
    scalar k2min = 3.0;
    scalar k2max = 10.0;

    if (3.0*x <= 1.0)
    {
        return k2min;
    }
    else if (3.0*x <= 2.0)
    {
        return k2min*sqr(2.0 - 3.0*x)*(6.0*x - 1.0)
            + k2max*sqr(3.0*x - 1.0)*(5.0 - 6.0*x);
    }

    return k2max;
}

int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("setMoments1DDict");
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

            forAll(moments, mi)
            {
                moments[mi][celli] = 0.0;
            }

            scalar alpha
                = 3.5 + 1.5*Foam::sin(2.0*Foam::constant::mathematical::pi*x);

            scalar beta
                = 3.5 - 1.5*Foam::cos(2.0*Foam::constant::mathematical::pi*x);

            moments[0][celli] = 1.0;

            for (label mi = 1; mi < nMoments; mi++)
            {
                moments[mi][celli] = moments[mi - 1][celli]
                    *(alpha + scalar(mi - 1))/(alpha + beta + scalar(mi - 1));
            }

            scalar coeffLeft
                = 0.5*
                (
                    1.0 + Foam::tanh(Foam::tan(Foam::constant::mathematical::pi
                    *(2.0*x - 0.5)))
                );

            scalar coeffRight =
                0.5*(1.0 + Foam::tanh(Foam::tan(Foam::constant::mathematical::pi
                *(-2.0*x + 1.5))));

            for (label mi = 0; mi < nMoments; mi++)
            {
                if (x <= 0.5)
                {
                    moments[mi][celli] *= coeffLeft;
                }
                else
                {
                    moments[mi][celli] *= coeffRight;
                }
            }
        }
    }
    else if (initType == "bimodal")
    {
        forAll(cellCentres, celli)
        {
            scalar x = cellCentres[celli].x();

            scalar w1, w2, w3, lambda2, k2, xi1;

            w1 = 0.0;
            w2 = 0.0;
            w3 = 0.0;

            if ((x >= 0.0) && (x <= 1.0))
            {
                w1 = 16.0*sqr(x)*sqr(1.0 - x);
            }


            if ((4.0*x >= 1.0) && (x <= 1.0))
            {
                w2 = 256.0*sqr(4.0*x - 1.0)*sqr(1.0 - x)/81.0;
            }

            w3 = weights2(x);
            lambda2 = parLambda2(x);
            k2 = park2(x);
            xi1 = 0.02;

            forAll(moments, mi)
            {
                moments[mi][celli]
                    = w1*pow(xi1, mi) + w2*pow(2.0*xi1, mi)
                      + w3*pow(lambda2, mi)*(std::tgamma(1.0 + scalar(mi)/k2));
            }
        }
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
