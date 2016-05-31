/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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
    Foam::populationBalanceModels::univariatePopulationBalance

Description
    Solve a univariate population balance equation with the extended quadrature
    method of moments.

SourceFiles
    univariatePopulationBalance.C

\*---------------------------------------------------------------------------*/

#ifndef univariatePopulationBalance_H
#define univariatePopulationBalance_H

#include "univariatePDFTransportModel.H"
#include "populationBalanceModel.H"
#include "aggregationKernel.H"
#include "breakupKernel.H"
#include "daughterDistribution.H"
#include "growthModel.H"
#include "diffusionModel.H"
#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{

/*---------------------------------------------------------------------------*\
                Class univariatePopulationBalance Declaration
\*---------------------------------------------------------------------------*/

class univariatePopulationBalance
:
    public univariatePDFTransportModel, public populationBalanceModel
{
    // Private data

        //- Name of the univariatePopulationBalance
        const word name_;

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

        autoPtr
        <
            Foam::populationBalanceSubModels::nucleationModel
        > nucleationModel_;

    // Private member functions

        //- Calculate aggregation source terms
        tmp<volScalarField> aggregationSource
        (
            const volUnivariateMoment& moment
        );

        //- Calculate Breakup source term
        tmp<volScalarField> breakupSource(const volUnivariateMoment& moment);

        //- Diffusion term
        virtual tmp<fvScalarMatrix> momentDiffusion
        (
            const volUnivariateMoment& moment
        );

        //- Convect moments in phase space
        virtual tmp<volScalarField> phaseSpaceConvection
        (
            const volUnivariateMoment& moment
        );

        //- Calculate source terms
        virtual tmp<fvScalarMatrix> momentSource
        (
            const volUnivariateMoment& moment
        );


public:

    //- Runtime type information
    TypeName("univariate");


    // Constructors

        //- Construct from components
        univariatePopulationBalance
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~univariatePopulationBalance();

    // Member Functions

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
