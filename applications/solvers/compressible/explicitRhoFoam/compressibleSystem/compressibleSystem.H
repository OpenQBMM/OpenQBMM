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
    Foam::compressibleSystem

Description
    Class used to handle encoding/decoding of primative and conservative
    variables used in compressible flows. Also includes funcitons to calculate
    the conservative fluxes

SourceFiles
    compressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleSystem_H
#define compressibleSystem_H

#include "IOdictionary.H"
#include "rhoThermo.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fluxIntegrator;
class fluxFunction;

/*---------------------------------------------------------------------------*\
                  Class compressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class compressibleSystem
{
protected:
    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Rho thermo class
        autoPtr<rhoThermo> thermoPtr_;

        //- Density field
        volScalarField rho_;

        //- velocity field
        volVectorField U_;

        //- Volumetric flux
        autoPtr<surfaceScalarField> phiPtr_;

        //- reference to thermodynamic pressure
        volScalarField& p_;

        //- Total internal energy
        volScalarField E_;

        //- Total enthalpy
        volScalarField H_;

        //  Conservative quantities

            //- Momentum
            volVectorField rhoU_;

            //- Total energy
            volScalarField rhoE_;

        //  Conservative fluxes

            //- Mass flux
            surfaceScalarField massFlux_;

            //- Momentum flux
            surfaceVectorField momentumFlux_;

            //- Energy flux
            surfaceScalarField energyFlux_;

        //- Ode steping data

            //- Flux integrator
            autoPtr<fluxIntegrator> integrator_;

            //- fields to store for ode time stepping
            labelList storedFieldIndexes_;
            labelList storedDeltaIndexes_;

            //- Changes in mass
            PtrList<volScalarField> rhos_;
            PtrList<volScalarField> deltaRhos_;

            //- Changes in momentum
            PtrList<volVectorField> rhoUs_;
            PtrList<volVectorField> deltaRhoUs_;

            //- Changes in energy
            PtrList<volScalarField> rhoEs_;
            PtrList<volScalarField> deltaRhoEs_;

        //- Flux function
        autoPtr<fluxFunction> fluxFunction_;


public:

    // Constructor
    compressibleSystem(const fvMesh& mesh);

    //- Destructor
    virtual ~compressibleSystem();


    // Member Functions

        //- Return speed of sound
        tmp<volScalarField> speedOfSound() const;

        //- Const reference to thermo
        inline const rhoThermo& thermo() const;

        //- Non-const reference to thermo
        inline rhoThermo& thermo();

        //- Const reference to density
        inline const volScalarField& rho() const;

        //- Non-const reference to density
        inline volScalarField& rho();

        //- Const reference to total internal energy
        inline const volScalarField& E() const;

        //- Non-const reference to total energy
        inline volScalarField& E();

        //- Const reference to total enthalpy
        inline const volScalarField& H() const;

        //- Non-const reference to total enthalpy
        inline volScalarField& H();

        //- Const reference to velocity
        inline const volVectorField& U() const;

        //- Non-const reference to velocity
        inline volVectorField& U();

        //- Const reference to volumetric flux
        inline const surfaceScalarField& phi() const;

        //- Non-const reference to volumetric Flux
        inline surfaceScalarField& phi();

        //- Const reference to thermodynamic pressure
        inline const volScalarField& p() const;

        //- Non-const reference to thermodynamic pressure
        inline volScalarField& p();

        //- Const reference to momentum
        inline const volVectorField& rhoU() const;

        //- Non-const reference to momentum
        inline volVectorField& rhoU();

        //- Const reference to total internal energy
        inline const volScalarField& rhoE() const;

        //- Non-const reference to total internal energy
        inline volScalarField& rhoE();

        //- Constant reference to mass flux
        inline const surfaceScalarField& massFlux() const;

        //- Non-constant reference to mass flux
        inline surfaceScalarField& massFlux();

        //- Constant reference to momentum flux
        inline const surfaceVectorField& momentumFlux() const;

        //- Non-constant reference to momentum flux
        inline surfaceVectorField& momentumFlux();

        //- Constant reference to energy flux
        inline const surfaceScalarField& energyFlux() const;

        //- Non-constant reference to energy flux
        inline surfaceScalarField& energyFlux();


    // Edit functions

        //- Set the number of steps used in ode solution to hyperbolic terms
        void setNSteps
        (
            const boolList& storeFields,
            const boolList& storeDeltas
        );

        //- Advect with conservative fluxes
        void advect
        (
            const label stepi,
            const scalarList& conservedVariablesCoeffs,
            const scalarList& fluxCoeffs,
            const dimensionedScalar& deltaT,
            const dimensionedVector& g
        );

        //- Transport conservative variables
        void integrateFluxes(const dimensionedVector& g);

        //- Update fluxes
        void updateFluxes();

        //- Correct primatives based on internal energy/enthalpy
        void correctThermo();

        //- Calculate primative variables
        void calcConservativeVariables();

        //- Calculate conservative variables
        void calcPrimitiveVariables();

        //- Read base phaseProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
