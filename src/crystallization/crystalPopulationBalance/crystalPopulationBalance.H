/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
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
    Foam::populationBalanceModels::crystalPopulationBalance

Description
    Solve a univariate population balance equation with the extended quadrature
    method of moments.

SourceFiles
    crystalPopulationBalance.C

    Considering the common number density unit of aerosol particle is  /cm^3,
    here we set the default unit for weight/m0 to  /cm^3. 
    This also reduces the number maganitude, which helps with the calculation . 
    However, in the nucleation and growth model,  necessary unit conversion is conducted. 
    The default size is 1nm. 

\*---------------------------------------------------------------------------*/

#ifndef crystalPopulationBalance_H
#define crystalPopulationBalance_H

#include "univariatePDFTransportModel.H"
#include "populationBalanceModel.H"

#include "nucleationModel.H"
#include "growthModel.H"
#include "aggregationKernel.H"
#include "breakupKernel.H"
#include "diffusionModel.H"
#include "depositionModel.H"
#include "solutionSaturationModel.H"
#include "realizableOdeSolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{

/*---------------------------------------------------------------------------*\
                Class crystalPopulationBalance Declaration
\*---------------------------------------------------------------------------*/

class crystalPopulationBalance
:
    public univariatePDFTransportModel,
    public populationBalanceModel,
    private realizableOdeSolver<volScalarMoment, volScalarNode>
{
    typedef realizableOdeSolver<volScalarMoment, volScalarNode> odeType;

    // Private dat

        Switch fixedSize_; 

        // minimum aerosol paricle number density, 10 /cm^3
        const scalar minM0_;

        //- Switch indicating the status of the nucleation model
        Switch nucleation_;

        //- Switch indicating the status of the growth model
        Switch growth_;

        //- Switch indicating the status of the aggregation model
        Switch aggregation_;

        //- Switch indicating the status of the breakup model
        Switch breakup_;

        //- Switch for deposition
        Switch deposition_;

        //- Switch for deposition
        Switch saturation_;

        //- Nucleation model
        autoPtr
        <
            Foam::populationBalanceSubModels::nucleationModel
        > nucleationModel_;

        //- Growth model function
        autoPtr
        <
            Foam::populationBalanceSubModels::growthModel
        > growthModel_;

        //- Diffusion model
        autoPtr
        <
            Foam::populationBalanceSubModels::diffusionModel
        > diffusionModel_;

        
        //- Aggregation kernel model
        autoPtr
        <
            Foam::populationBalanceSubModels::aggregationKernel
        > aggregationKernel_;

        //- Breakup kernel model
        autoPtr
        <
            Foam::populationBalanceSubModels::breakupKernel
        > breakupKernel_;

        //- deposition kernel
        autoPtr
        <
            populationBalanceSubModels::depositionModel
        > depositionModel_;


        //- saturationModel
        autoPtr
        <
            solutionSaturationModel
        > saturationModel_;

        //- Crystal density (kg/m3)
        dimensionedScalar rhop_;
        
        //- Shape factor for crystal (sphere: pi/6)
        dimensionedScalar shapeFactor_;

        //- Prandtl number
        dimensionedScalar Pr_;
        
        //- Turbulent Prandtl number
        dimensionedScalar Prt_;
        
        //- Schmidt number
        dimensionedScalar Sct_;
        
        //- Molar diffusivity
        dimensionedScalar molarDiff_;
        
        //- Number mean crystal length L10 = m1/m0
        volScalarField L10_;

        //- Sauter mean crystal length L32 = m3/m2
        volScalarField L32_;

        //- Legacy mean crystal size output. Prefer L10_ or L32_.
        volScalarField d_;

        dimensionedScalar d_nucleation_;

    // Private member functions

        //- Sub Models preUpdate
        virtual void subModelsPreUpdate();

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volScalarMoment& moment
        );

        //- Calculate explicit source terms with realizable ODE solver
        virtual void explicitMomentSource();

        //- Are moments updated due to sources
        virtual bool solveMomentSources() const;

        //- Is the realizable ode solver used to solve sources
        virtual bool solveMomentOde() const;

        //- Update explicit collison source terms
        virtual void updateCellMomentSource(const label celli);

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& nodes,
            const label environment
        );

        virtual void checkCorrectMoments(); 

        void updateCharacteristicLengths();

        void calcSpeciesTransfer();

public:

    //- Runtime type information
    TypeName("crystal");


    // Constructors

        //- Construct from components
        crystalPopulationBalance
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~crystalPopulationBalance();

    // Member Functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const;

        //- Return the courant number based on abscissae
        virtual scalar CoNum() const;

        //- Solve population balance equation
        void solve();

        //- Read
        virtual bool readIfModified();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace univariatePDFTransportModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
