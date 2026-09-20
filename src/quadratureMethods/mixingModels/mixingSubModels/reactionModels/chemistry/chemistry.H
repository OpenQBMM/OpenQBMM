/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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
    Foam::mixingSubModels::reactionModels::chemistry

Description
    The chemistry of OpenFOAM solved in every environment of the mixture
    fraction.

    The compositions conditioned on the mixture fraction are the mass
    fractions of the species of the reacting thermo of the case, in the
    order of its species list: the composition j of the quadrature is the
    species j, and the quadrature must carry as many compositions as the
    thermo has species, which the model checks. The moment fields keep
    their index names; the species are named only where OpenFOAM names
    them, in the thermophysical properties, the reactions and the thermo
    files of the case, which are the same files a reacting solver reads.

    An environment reacts as one cell of a reacting solver does, on its own
    state: the concentrations of its species, c_j = rho_i Y_ij / W_j, with
    rho_i the density of the environment at its own temperature and
    composition and the pressure of the cell, its temperature T_i and the
    pressure. The state is advanced over the time step of the case by the
    ODE solver of the odeCoeffs sub-dictionary of chemistryProperties,
    working on the derivatives the chemistry model of OpenFOAM gives for
    its species and its temperature: the reaction rates of the mechanism,
    and the temperature at constant pressure with no heat exchanged, so
    that the enthalpy of an environment is what the reaction conserves.
    Which chemistry model OpenFOAM builds, and with which chemistry solver,
    is what chemistryType selects, as in any reacting case; the
    integration here always uses odeCoeffs, whatever that solver is, and
    keeps for every environment of every cell its own estimate of the
    chemical time step, started at initialChemicalTimeStep and never above
    maxChemicalTimeStep. An environment below Treact does not react. The
    mass fractions written back are the concentrations normalised by the
    mass they carry, so that their sum stays one whatever round-off the
    integration leaves.

    The temperature of an environment is a closure of the model,
    selected by temperature:
    - isothermal: every environment is at the temperature of the cell,
      which the model leaves as the case gives it, and stays at it while
      it reacts: the temperature equation of the chemistry model is
      switched off, so that the heat of reaction changes nothing;
    - linearEnthalpy: the enthalpy of an environment is linear in the
      mixture fraction between the enthalpies of the streams at xi = 0 and
      xi = 1, which is what mixing alone gives when nothing else changes
      the enthalpy: adiabatic walls, no radiation, and heat and species
      diffusing alike. The streams are given by their temperature and
      composition, and the temperature of an environment is the one at
      which its composition has that enthalpy. The absolute enthalpy is
      used, so that the heat of reaction is what raises the temperature.
    A closure carrying the enthalpy of every environment as a coordinate
    of the distribution, mixed with the compositions and untouched by
    reaction, is not here yet.

    The temperature of every environment is kept in a field named after
    the node, T.node0.mixing for the first node of the distribution
    mixing, beside the weight and the abscissae of the node, which hold
    its mixture fraction and its compositions: weight.node0.mixing,
    abscissa0.node0.mixing for the mixture fraction, abscissa1.node0.mixing
    for the first composition. None of them is written with the time
    directories, but all are registered, so the probes function object
    samples them by name, and the writeObjects function object writes
    them when a case wants them on disk.

    Once the environments of a cell have reacted, the moments are rebuilt
    from the nodes, exactly, and the mean state of the cell is written to
    the thermo, so that the properties of the fluid follow the quadrature:
    the mass fraction fields of the species receive the mean composition,
    the mean over the environments weighted by the nodes, and, with
    linearEnthalpy, the energy field receives the enthalpy of the mean
    temperature, the mean of the temperatures of the environments. The
    solver updates the thermo after the mixing model has been solved.
    The mean density of a variable-density flow, the mean over the
    environments of their densities, is not written here.

    The cost is that of one cell of a reacting solver for every
    environment, so a cell costs n times what it costs reactingFoam with
    the same mechanism and the same ODE solver.

    The case needs what a reacting solver needs: a thermo with a
    reactingMixture, whose species fields are initialised from Ydefault
    or one file per species, the reactions and thermo files it names, and
    chemistryProperties. The model itself reads:

    \verbatim
    reaction
    {
        reactionModel   chemistry;
        temperature     linearEnthalpy;     // isothermal | linearEnthalpy
        streams                             // with linearEnthalpy
        {
            xi0 { T 293; Y { O2 0.23; N2 0.77; } }
            xi1 { T 293; Y { CH4 1; } }
        }
    }
    \endverbatim

    The mass fractions of a stream are given by species name, those not
    given are zero, and they must sum to one.

    With the debug switch chemistry on, the state of every environment is
    printed before and after it reacts.

SourceFiles
    chemistry.C

\*---------------------------------------------------------------------------*/

#ifndef chemistry_H
#define chemistry_H

#include "reactionModel.H"
#include "basicThermo.H"
#include "basicSpecieMixture.H"
#include "basicChemistryModel.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                          Class chemistry Declaration
\*---------------------------------------------------------------------------*/

class chemistry
:
    public reactionModel
{
public:

    //- How the temperature of an environment is found
    enum class temperatureClosure
    {
        isothermal,
        linearEnthalpy
    };

    //- Names of the temperature closures
    static const Enum<temperatureClosure> temperatureClosureNames;


private:

    // Private classes

        //- The ODE system of the chemistry model with its temperature
        //  held: the state of an environment that reacts at the
        //  temperature it is given
        class isothermalSystem
        :
            public ODESystem
        {
            //- The system of the chemistry model
            const ODESystem& system_;

            //- Index of the temperature in the state
            const label TIndex_;

        public:

            //- Construct from the system and the index of its temperature
            isothermalSystem(const ODESystem& system, const label TIndex)
            :
                system_(system),
                TIndex_(TIndex)
            {}

            //- Number of equations
            virtual label nEqns() const
            {
                return system_.nEqns();
            }

            //- The derivatives of the system, none for the temperature
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- The Jacobian of the system, the temperature depending on
            //  nothing
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;
        };


    // Private data

        //- The thermo of the case, which owns the species
        basicThermo& thermo_;

        //- The composition of the thermo: species, their data and fields
        basicSpecieMixture& composition_;

        //- The chemistry model of OpenFOAM, built on the thermo
        autoPtr<basicChemistryModel> chemistry_;

        //- The chemistry model as the ODE system of the state of an
        //  environment: its concentrations, temperature and pressure
        const ODESystem& system_;

        //- The temperature closure
        const temperatureClosure temperature_;

        //- That system with its temperature held, for isothermal
        //  environments
        autoPtr<isothermalSystem> isothermalSystem_;

        //- The ODE solver of odeCoeffs, integrating the system the closure
        //  calls for
        autoPtr<ODESolver> odeSolver_;

        //- Temperature of the stream at xi = 0
        scalar T0_;

        //- Composition of the stream at xi = 0
        scalarField Y0_;

        //- Temperature of the stream at xi = 1
        scalar T1_;

        //- Composition of the stream at xi = 1
        scalarField Y1_;

        //- Temperature below which an environment does not react
        const scalar Treact_;

        //- Initial estimate of the chemical time step
        const scalar deltaTChemIni_;

        //- Largest chemical time step
        const scalar deltaTChemMax_;

        //- Estimate of the chemical time step of every environment of
        //  every cell
        List<scalarField> subDeltaT_;

        //- The temperature of every environment
        PtrList<volScalarField> T_;

        //- The state of an environment: concentrations, temperature and
        //  pressure
        scalarField cTp_;

        //- The composition of an environment
        scalarField Y_;


    // Private member functions

        //- The composition of a thermo that has species
        static basicSpecieMixture& compositionOf(basicThermo& thermo);

        //- The chemistry model of OpenFOAM for the thermo
        static autoPtr<basicChemistryModel> newChemistry(basicThermo& thermo);

        //- The chemistry model as an ODE system
        static const ODESystem& systemOf(const basicChemistryModel& model);

        //- Read a stream: its temperature and its composition by species
        void readStream
        (
            const dictionary& dict,
            scalar& T,
            scalarField& Y
        ) const;

        //- Absolute enthalpy of a composition [J/kg]
        scalar enthalpy
        (
            const scalarField& Y,
            const scalar p,
            const scalar T
        ) const;

        //- Heat capacity at constant pressure of a composition [J/kg/K]
        scalar Cp
        (
            const scalarField& Y,
            const scalar p,
            const scalar T
        ) const;

        //- Density of a composition [kg/m3]
        scalar density
        (
            const scalarField& Y,
            const scalar p,
            const scalar T
        ) const;

        //- Temperature of an environment by the closure
        scalar nodeTemperature
        (
            const scalar xi,
            const scalarField& Y,
            const scalar Tcell,
            const scalar p
        ) const;


public:

    //- Runtime type information
    TypeName("chemistry");


    // Constructors

        //- Construct from components
        chemistry
        (
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            const label nCompositions
        );


    //- Destructor
    virtual ~chemistry();


    // Member Functions

        //- The temperature closure
        temperatureClosure temperature() const
        {
            return temperature_;
        }

        //- The temperature of every environment
        const PtrList<volScalarField>& T() const
        {
            return T_;
        }

        //- Whether the chemistry of the case is on
        virtual bool reacts() const;

        //- React every environment of every cell over the time step of
        //  the case, rebuild the moments and write the mean state of the
        //  cell to the thermo
        virtual void react(velocityQuadratureApproximation& quadrature);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels
} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
