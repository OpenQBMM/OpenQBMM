/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
     \\/     M anipulation  |
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
    Foam::populationBalanceModels::mixingPopulationBalance

Description
    Solve a coupled mixing-population balance problem accounting for multiple
    environments.

SourceFiles
    mixingPopulationBalance.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPopulationBalance_H
#define mixingPopulationBalance_H

#include "univariatePDFTransportModel.H"
#include "populationBalanceModel.H"

#include "turbulentMixing.H"

#include "aggregationKernel.H"
#include "breakupKernel.H"
#include "daughterDistribution.H"
#include "growthModel.H"
#include "diffusionModel.H"
#include "nucleationModel.H"
#include "environmentMixingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{

/*---------------------------------------------------------------------------*\
                    Class mixingPopulationBalance Declaration
\*---------------------------------------------------------------------------*/

class mixingPopulationBalance
:
    public populationBalanceModel
{
    // Private data

        //- Absolute tolerance for realizable ODE solver
        scalar ATol_;

        //- Relative tolerance for realizable ODE solver
        scalar RTol_;

        //- Time-step change factor
        scalar fac_;

        //- Minimum time-step change factor
        scalar facMin_;

        //- Maximum time-step change factor
        scalar facMax_;

        //- Minimum local step in realizable ODE solver
        label minLocalDt_;

        //- Name of the mixingPopulationBalance
        const word name_;

        //- Turbulent mixing model
        autoPtr
        <
            Foam::PDFTransportModels::mixingModels::turbulentMixing
        > mixingModel_;

        //- Reference to field of weight of environment one
        const volScalarField& p1_;

        //- Reference to field of abscissa of environment one
        const volScalarField& xi1_;

        //- Reference to field of weight of environment two
        const volScalarField& p2_;

        //- Reference to field of abscissa of environment two
        const volScalarField& xi2_;

        //- Quadrature for the mean moments
        univariateQuadratureApproximation meanMomentsQuadrature_;

        //- Quadrature for the mean moment variance
        univariateQuadratureApproximation meanMomentsVarianceQuadrature_;

        //- Moment advection for mean moments
        autoPtr<univariateMomentAdvection> meanMomentsAdvection_;

        //- Moment advection for variance of mean moments
        autoPtr<univariateMomentAdvection> meanMomentsVarianceAdvection_;

        //- Quadrature approximation of the PBE in the first environment
        univariateQuadratureApproximation envOneQuadrature_;

        //- Quadrature approximation of the PBE in the second environment
        univariateQuadratureApproximation envTwoQuadrature_;

        //- Reference to moments in environment one
        volUnivariateMomentFieldSet& mEnvOne_;

        //- Reference to moments in environment two
        volUnivariateMomentFieldSet& mEnvTwo_;

        //- Switch indicating the status of the aggregation model
        Switch aggregation_;

        //- Switch indicating the status of the breakup model
        Switch breakup_;

        //- Switch indicating the status of the growth model
        Switch growth_;

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

        //- Daughter distribution function
        autoPtr
        <
            Foam::populationBalanceSubModels::daughterDistribution
        > daughterDistribution_;

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

        //- Nucleation model
        autoPtr
        <
            Foam::populationBalanceSubModels::nucleationModel
        > nucleationModel_;

        //- Environment mixing model
        autoPtr
        <
            Foam::populationBalanceSubModels::environmentMixingModel
        > envMixingModel_;


    // Private member functions

        //- Calculate explicit source terms with realizable ODE solver
        void explicitMomentSource
        (
            univariateQuadratureApproximation& quadrature
        );

        //- Calculate moments of the PBE in each environment
        void calcEnvironmentMoments();

        //- Calculated mixed moments
        void calcMixedMoments();

        //- Calculate aggregation source terms
        scalar aggregationSource
        (
            const label momentOrder,
            const label celli,
            const PtrList<volScalarNode>& nodes
        );

        //- Calculate Breakup source term
        scalar breakupSource
        (
            const label momentOrder,
            const label celli,
            const PtrList<volScalarNode>& nodes
        );

        //- Diffusion term
        virtual tmp<fvScalarMatrix> momentDiffusion
        (
            const volUnivariateMoment& moment
        );

        //- Convect moments in phase space
        virtual scalar phaseSpaceConvection
        (
            const label momentOrder,
            const label celli,
            const PtrList<volScalarNode>& nodes
        );

        //- Calculate source term for moment equation in the specified cell
        scalar cellMomentSource
        (
            label momentOrder,
            label celli,
            const PtrList<volScalarNode>& nodes
        );


public:

    //- Runtime type information
    TypeName("mixingPbe");


    // Constructors

        //- Construct from components
        mixingPopulationBalance
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~mixingPopulationBalance();

    // Member Functions

        //- Return the maximum Courant number ensuring moment realizability
        scalar realizableCo();

        //- Solve population balance equation
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace univariatePDFTransportModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
