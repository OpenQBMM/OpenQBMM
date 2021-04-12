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
    Foam::AUSM

Description
    Computes fluxes using AUSM scheme

SourceFiles
    AUSMFlux.C

\*---------------------------------------------------------------------------*/

#ifndef AUSMFlux_H
#define AUSMFlux_H

#include "fluxFunction.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxFunctions
{

/*---------------------------------------------------------------------------*\
                  Class fluxFunction Declaration
\*---------------------------------------------------------------------------*/

class AUSM : public fluxFunction
{
public:

    //- Runtime type information
    TypeName("AUSM");

    // Costructor
    AUSM(const fvMesh& mesh);


    //- Destructor
    virtual ~AUSM();


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
            const volScalarField& a
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace fluxFunctions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
