/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 Alberto Passalacqua
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

Class
    Foam::fluxFunction

Description
    Interpolates conservative and primative variables and updates fluxes

SourceFiles
    fluxFunction.C
    newFluxFunction.C

\*---------------------------------------------------------------------------*/

#ifndef fluxFunction_H
#define fluxFunction_H

#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fluxFunction Declaration
\*---------------------------------------------------------------------------*/

class fluxFunction
{
protected:
    // Protected Data

        //- Refrence to compressibleSinglePhaseSystem
        const fvMesh& mesh_;

        const dictionary& dict_;

        //- Owner interpolation
        surfaceScalarField own_;

        //- Neighbour interpolation
        surfaceScalarField nei_;

        //- Residual alpha
        dimensionedScalar residualAlpha_;

        //- Cutoff Mach number
        dimensionedScalar cutoffMa_;


    //- Protected functions

        //- Return interpolation scheme name
        static word schemeName(const word& fieldName)
        {
            return "reconstruct(" + fieldName + ")";
        }

public:

    //- Runtime type information
    TypeName("fluxFunction");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxFunction,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructor
    fluxFunction(const fvMesh& mesh);

    // Selector
    static autoPtr<fluxFunction> New(const fvMesh& mesh);

    //- Destructor
    virtual ~fluxFunction();


    // Public Functions

        //- Update fluxes
        virtual void updateFluxes
        (
            surfaceScalarField& massFlux,
            surfaceVectorField& momentumFlux,
            surfaceScalarField& energyFlux,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& H,
            const volScalarField& p,
            const volScalarField& gamma
        ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
