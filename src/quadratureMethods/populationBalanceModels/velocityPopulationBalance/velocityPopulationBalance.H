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
    Copyright (C) 2019 Alberto Passalacqua
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

#include "velocityPDFTransportModel.H"
#include "populationBalanceModel.H"
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
    public velocityPDFTransportModel,
    public populationBalanceModel,
    public realizableOdeSolver<volVelocityMoment, volVelocityNode>
{
protected:

    typedef realizableOdeSolver<volVelocityMoment, volVelocityNode> odeType;

    // protected data

        //- Switch indicating the status of the aggregation model
        Switch collision_;

        //- Aggregation kernel model
        autoPtr
        <
            Foam::populationBalanceSubModels::collisionKernel
        > collisionKernel_;


    // Private member functions

        //- Update implicit collison source terms
        virtual void updateImplicitMomentSource();

         //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volVelocityMoment& moment
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
            const velocityQuadratureApproximation& quadrature,
            const label enviroment = 0
        );


public:

    //- Runtime type information
    TypeName("velocity");


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

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const;

        //- Return the maximum Courant number based on velocities
        virtual scalar CoNum() const;

        //- Solve population balance equation
        virtual void solve();

        //- Read
        virtual bool readIfModified();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace velocityPDFTransportModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
