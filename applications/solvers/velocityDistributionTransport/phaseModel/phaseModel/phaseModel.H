/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2019-04-29 Jeff Heylmun:    Simplified model
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

Class
    Foam::phaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    phaseModel.C
    newPhaseModel.C
    phaseModels.C

\*---------------------------------------------------------------------------*/

#ifndef pbePhaseModel_H
#define pbePhaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "transportModel.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Phase>
class PhaseCompressibleTurbulenceModel;


/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField,
    public transportModel
{
protected:
    // Protected data

        //- Name of phase
        word name_;

        //- Phase dictionary
        dictionary phaseDict_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Optional maximum phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        //- Thermophysical properties
        autoPtr<rhoThermo> thermo_;

        //- Mean velocity
        volVectorField U_;

        //- Diameter
        volScalarField d_;

        //- Volumetric flux of the phase
        autoPtr<surfaceScalarField> phiPtr_;

        //- Volume fraction flux
        surfaceScalarField alphaPhi_;

        //- Mass flux
        surfaceScalarField alphaRhoPhi_;

        //- Turbulence model
        autoPtr<PhaseCompressibleTurbulenceModel<phaseModel>> turbulence_;


public:

    //- Runtime type information
    ClassName("phaseModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseModel,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& phaseName
            ),
            (mesh, dict, phaseName)
        );

    // Constructors
        phaseModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName = word::null
        );

         //- Return clone
        autoPtr<phaseModel> clone() const;


    // Selectors

        static autoPtr<phaseModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the name of this phase
        const word& name() const
        {
            return name_;
        }

        //- Return the number of nodes
        virtual label nNodes() const
        {
            return 1;
        }

        //- Return the turbulence model
        const PhaseCompressibleTurbulenceModel<phaseModel>& turbulence() const;

        //- Return non-const access to the turbulence model
        //  for correction
        PhaseCompressibleTurbulenceModel<phaseModel>& turbulence();

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        const dimensionedScalar& residualAlpha() const
        {
            return residualAlpha_;
        }

        //- Constant access to alpha field for nodei
        virtual tmp<volScalarField> volumeFraction(const label nodei = -1) const
        {
            return *this;
        }

        //- Return diameter of nodei
        virtual tmp<volScalarField> d(const label nodei = -1) const
        {
            return d_;
        }

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return the thermophysical model
        const rhoThermo& thermo() const
        {
            return thermo_();
        }

        //- Return non-const access to the thermophysical model
        //  for correction
        rhoThermo& thermo()
        {
            return thermo_();
        }

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return thermo_->nu();
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return thermo_->nu(patchi);
        }

        //- Return the laminar dynamic viscosity
        tmp<volScalarField> mu() const
        {
            return thermo_->mu();
        }

        //- Return the laminar dynamic viscosity for patch
        tmp<scalarField> mu(const label patchi) const
        {
            return thermo_->mu(patchi);
        }

        //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
        tmp<volScalarField> alpha() const
        {
            return thermo_->alpha();
        }

        //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
        tmp<scalarField> alpha(const label patchi) const
        {
            return thermo_->alpha(patchi);
        }

        //- Thermal diffusivity for temperature of mixture [J/m/s/K]
        tmp<scalarField> kappa(const label patchi) const
        {
            return thermo_->kappa(patchi);
        }

        //- Thermal diffusivity for temperature of mixture
        //  for patch [J/m/s/K]
        tmp<volScalarField> kappa() const
        {
            return thermo_->kappa();
        }

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        tmp<volScalarField> alphahe() const
        {
            return thermo_->alphahe();
        }

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        tmp<scalarField> alphahe(const label patchi) const
        {
            return thermo_->alphahe(patchi);
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [J/m/s/K]
        tmp<volScalarField> kappaEff
        (
            const volScalarField& alphat
        ) const
        {
            return thermo_->kappaEff(alphat);
        }

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [J/m/s/K]
        tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const
        {
            return thermo_->kappaEff(alphat, patchi);
        }

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const
        {
            return thermo_->alphaEff(alphat);
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const
        {
            return thermo_->alphaEff(alphat, patchi);
        }

        //- Return the specific heat capacity
        tmp<volScalarField> Cp() const
        {
            return thermo_->Cp();
        }

        //- Return the density
        const volScalarField& rho() const
        {
            return thermo_->rho();
        }

        //- Constant access to the mean velocity
        virtual const volVectorField& U(const label nodei = -1) const
        {
            return U_;
        }

        //- Non-const access to the mean velocity
        virtual volVectorField& U(const label nodei = -1)
        {
            return U_;
        }

        //- Deviation of velocity nodei from the mean
        virtual tmp<volVectorField> Vs(const label nodei = -1) const;

        //- Constant access to the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        //- Non-const access to the volumetric flux
        surfaceScalarField& phi()
        {
            return phiPtr_();
        }

        //- Constant access to the volume fraction flux
        const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }

        //- Non-const access to the volume fraction flux
        surfaceScalarField& alphaPhi()
        {
            return alphaPhi_;
        }

        //- Constant access to the mass flux
        const surfaceScalarField& alphaRhoPhi() const
        {
            return alphaRhoPhi_;
        }

        //- Non-const access to the mass flux
        surfaceScalarField& alphaRhoPhi()
        {
            return alphaRhoPhi_;
        }

        //- Solve the underlying population balance model
        virtual void solve()
        {}

        //- Correct
        virtual void correct()
        {}

        //- Relative transport of moments, not used in standard phase model
        virtual void relativeTransport()
        {
            return;
        }

        //- Average transport of moments, not used in standard phase model
        virtual void averageTransport()
        {
            return;
        }

        //- Solve change in size moments based on breakup and coalesence
        virtual void solveSource()
        {
            return;
        }

        //- Return the maximum realizable courant number
        virtual scalar realizableCo() const
        {
            return 1.0;
        }

        //- Return the max courant number
        virtual scalar CoNum() const
        {
            return 0.0;
        }

        //- Read phase properties dictionary
        virtual bool read()
        {
            return true;
        }

        //- Read phase properties dictionary
        virtual bool read(const bool readOK)
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
