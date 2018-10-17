/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::kineticTheoryModel

Description
    Kinetic theory particle phase RAS model

    Reference:
    \verbatim
        "Derivation, implementation, and validation of computer simulation
         models for gas-solid fluidized beds",
        van Wachem, B.G.M.,
        Ph.D. Thesis, Delft University of Technology, Amsterdam, 2000.
    \endverbatim

    There are no default model coefficients.

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "phaseModel.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "radialModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"
#include "fvc.H"

#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
:
    public regIOobject
{
protected:

    // Protected data

            // Refrence to phase model
            const phaseModel& phase_;


        // Sub-models

            //- Run-time selected viscosity model
            autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;

            //- Run-time selected conductivity model
            autoPtr<kineticTheoryModels::conductivityModel> conductivityModel_;

            //- Run-time selected radial distribution model
            autoPtr<kineticTheoryModels::radialModel> radialModel_;

            //- Run-time selected granular pressure model
            autoPtr<kineticTheoryModels::granularPressureModel>
                granularPressureModel_;

            //- Run-time selected frictional stress model
            autoPtr<kineticTheoryModels::frictionalStressModel>
                frictionalStressModel_;


        // Kinetic Theory Model coefficients

            //- Use equilibrium approximation: generation == dissipation
            Switch equilibrium_;

            //- Coefficient of restitution
            dimensionedScalar e_;

            //- Maximum packing phase-fraction
            dimensionedScalar alphaMax_;

            //- Min value for which the frictional stresses are zero
            dimensionedScalar alphaMinFriction_;

            //- Residual phase fraction
            dimensionedScalar residualAlpha_;

            //- Maximum turbulent viscosity
            dimensionedScalar maxNut_;


        // Kinetic Theory Model Fields

            //- The granular energy/temperature
            volScalarField Theta_;

            //- The granular bulk viscosity
            volScalarField lambda_;

            //- The granular radial distribution
            volScalarField g0_;

            //- The granular "thermal" conductivity
            volScalarField kappa_;

            //- The particle viscosity
            volScalarField nu_;

            //- The frictional viscosity
            volScalarField nuFric_;


public:

    //- Runtime type information
    TypeName("kineticTheoryModel");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        kineticTheoryModel,
        dictionary,
        (
            const dictionary& dict,
            const phaseModel& phase
        ),
        (dict, phase)
    );


    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    // Selectors

        //- Select null constructed
        static autoPtr<kineticTheoryModel> New
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~kineticTheoryModel();


    // Member Functions

        //- Return the granular temperature
        virtual const volScalarField& Theta() const
        {
            return Theta_;
        }

        //- Return granular stress tensor
        virtual tmp<volSymmTensorField> Sigma() const;

        //- Return the granular temperature on a patch
        virtual const scalarField& ThetaPatch(const label patchi) const
        {
            return Theta_.boundaryField()[patchi];
        }

        //- Return the granular viscosity
        virtual const volScalarField& nu() const
        {
            return nu_;
        }

        //- Return the granular viscosity on a patch
        virtual const scalarField& nuPatch(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the effective viscosity (granular + frictional)
        virtual tmp<volScalarField> nuEff() const;

        //- Return the granular bulk viscosity
        virtual const volScalarField& lambda() const
        {
            return lambda_;
        }

        //- Return the granular "thermal" conductivity
        virtual const volScalarField& kappa() const
        {
            return kappa_;
        }

        //- Return the granular "thermal" conductivity on a patch
        virtual const scalarField& kappaPatch(const label patchi) const
        {
            return kappa_.boundaryField()[patchi];
        }

        //- Return the radial distribution
        virtual const volScalarField& g0() const
        {
            return g0_;
        }

        //- Return the radial distribution on a patch
        virtual const scalarField& g0Patch(const label patchi) const
        {
            return g0_.boundaryField()[patchi];
        }

        //- Return the radial distribution derivative w.r.t. the particulate
        //  particulate phase fraction
        virtual tmp<volScalarField> g0Prime() const;

        //- Return the phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<volScalarField> pPrime() const;

        //- Return the face-phase-pressure'
        // (derivative of phase-pressure w.r.t. phase-fraction)
        virtual tmp<surfaceScalarField> pPrimef() const;

        //- Return the maximum value of the phase fraction
        virtual const dimensionedScalar& alphaMax() const
        {
            return alphaMax_;
        }


    // Update functions

        //- Update kinetic theory quantities (turbulent correct)
        virtual void solve
        (
            const volScalarField& beta,
            const volScalarField& alpha,
            const volTensorField& gradU,
            const volSymmTensorField D
        ) = 0;

        //- Transport moments
        virtual void transportMoments()
        {}

        //- Max time-step
        virtual scalar maxUxDx() const
        {
            return 0.0;
        }

        //- Return fraction to be used for hydrodynamic flux
        virtual tmp<volScalarField> h2() const
        {
             return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "one",
                        phase_.mesh().time().timeName(),
                        phase_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    phase_.mesh(),
                    dimensionedScalar("1", dimless, 1.0)
                )
            );
        }

        //- Return interpolated fraction to be used for hydrodynamic flux
        virtual tmp<surfaceScalarField> h2f() const
        {
             return tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "one",
                        phase_.mesh().time().timeName(),
                        phase_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    phase_.mesh(),
                    dimensionedScalar("1", dimless, 1.0)
                )
            );
        }

        virtual tmp<volScalarField> ddtAlphaDilute() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "zero",
                        phase_.mesh().time().timeName(),
                        phase_.mesh(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    phase_.mesh(),
                    dimensionedScalar("0", inv(dimTime), 0.0)
                )
            );
        }


        //- Update viscosity and conductivity
        void update();


    // IO functions

        //- Read
        virtual bool read()
        {
            return true;
        }

        //- Write
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
