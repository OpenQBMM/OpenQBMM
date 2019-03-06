/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

        //- Nucleation model
        autoPtr
        <
            Foam::populationBalanceSubModels::nucleationModel
        > nucleationModel_;

    // Private member functions

        //- Calculate aggregation source terms
        scalar aggregationSource
        (
            const label& momentOrder,
            const label& celli
        );

        //- Calculate Breakup source term
        scalar breakupSource
        (
            const label& momentOrder,
            const label& celli
        );

        //- Diffusion term
        virtual tmp<fvScalarMatrix> momentDiffusion
        (
            const volMoment& moment
        );

        //- Convect moments in phase space
        virtual scalar phaseSpaceConvection
        (
            const label& momentOrder,
            const label& celli
        );

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volMoment& moment
        );

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            label& momentOrder,
            label& celli
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
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~univariatePopulationBalance();

    // Member Functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const;

        //- Return the courant number based on abscissae
        virtual scalar CoNum() const;

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
