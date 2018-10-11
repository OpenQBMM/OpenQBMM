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

#include "velocityPDFTransportModel.H"
#include "populationBalanceModel.H"
#include "collisionKernel.H"

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
    public velocityPDFTransportModel, public populationBalanceModel
{
    // Private data

        //- Name of the velocityPopulationBalance
        const word name_;

        //- Switch indicating the status of the aggregation model
        Switch collision_;

        //- Aggregation kernel model
        autoPtr
        <
            Foam::populationBalanceSubModels::collisionKernel
        > collisionKernel_;


    // Private member functions

        //- Return if collisions are on
        virtual bool collision() const;

        //- Update implicit collison source terms
        virtual void updateImplicitCollisionSource();

        //- Update explicit collison source terms
        virtual void updateExplicitCollisionSource(const label celli);

         //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitCollisionSource
        (
            const volVectorMoment& moment
        );

        virtual scalar explicitCollisionSource
        (
            const label momenti,
            const label celli
        );

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const label momentOrder,
            const label celli
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

        //- Solve mean tranport
        virtual void meanTransport
        (
            const surfaceScalarField& phi,
            const bool wallCollisions
        );

        //- Solve relative transport
        virtual void relativeTransport
        (
            const mappedPtrList<volVectorField>& Vs,
            const bool wallCollisions
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace velocityPDFTransportModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
