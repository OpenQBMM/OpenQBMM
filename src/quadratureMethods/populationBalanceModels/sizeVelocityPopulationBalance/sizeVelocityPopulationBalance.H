/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 Alberto Passalacqua
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
    sizeVelocityPopulationBalance.C

\*---------------------------------------------------------------------------*/

#ifndef sizeVelocityPopulationBalance_H
#define sizeVelocityPopulationBalance_H

#include "velocityPopulationBalance.H"
#include "collisionKernel.H"
#include "aggregationKernel.H"
#include "breakupKernel.H"
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
                Class velocityPopulationBalance Declaration
\*---------------------------------------------------------------------------*/

class sizeVelocityPopulationBalance
:
    public velocityPopulationBalance
{
    // Private data

        //- Switch indicating the status of the aggregation model
        Switch aggregation_;

        //- Switch indicating the status of the breakup model
        Switch breakup_;

        //- Switch indicating the status of the growth model
        Switch growth_;

        //- Switch indicating the status of the nucleation model
        Switch nucleation_;

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

//         //- Update explicit collison source terms
//         virtual void updateCellMomentSource(const label celli);

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volVelocityMoment& moment
        );

        //- Calculate explicit source terms with realizable ODE solver
        virtual void explicitMomentSource();

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature,
            const label enviroment = 0
        );


public:

    //- Runtime type information
    TypeName("sizeVelocity");


    // Constructors

        //- Construct from components
        sizeVelocityPopulationBalance
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~sizeVelocityPopulationBalance();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace velocityPDFTransportModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
