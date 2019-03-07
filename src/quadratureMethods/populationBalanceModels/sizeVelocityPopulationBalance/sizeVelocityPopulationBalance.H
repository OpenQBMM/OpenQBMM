/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Alberto Passalacqua
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
    Foam::populationBalanceModels::velocityPopulationBalance

Description
    Solve a velocity population balance equation with the hyperbolic
    conditional quadrature method of moments with collisions.

SourceFiles
    velocityPopulationBalance.C

\*---------------------------------------------------------------------------*/

#ifndef velocityPopulationBalance_H
#define velocityPopulationBalance_H

#include "velocityPopulationBalance.H"
#include "collisionKernel.H"
#include "realizableOdeSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{

/*---------------------------------------------------------------------------*\
                Class velocityPopulationBalance Declaration
\*---------------------------------------------------------------------------*/

class velocityPopulationBalance
:
    public velocityPopulationBalanceModel
{
    // Private data

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

//         //- Update implicit collison source terms
//         virtual void updateImplicitMomentSource();
//
//          //- Calculate implicit source terms
//         virtual tmp<fvScalarMatrix> implicitMomentSource
//         (
//             const volVelocityMoment& moment
//         );

        //- Calculate explicit source terms with realizable ODE solver
        virtual void explicitMomentSource();

//         //- Update explicit collison source terms
//         virtual void updateCellMomentSource(const label celli);

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const label momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature,
            const label enviroment = 0
        );


public:

    //- Runtime type information
    TypeName("sizeVelocity");


    // Constructors

        //- Construct from components
        velocityPopulationBalance
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~velocityPopulationBalance();

    // Member Functions

        //- Solve population balance equation
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace velocityPDFTransportModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
