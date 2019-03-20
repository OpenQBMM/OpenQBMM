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
    Foam::populationBalanceSubModels::collisionKernels::esBGKCollision

Description
    Extension of the BGK collisional model that includes the coefficient
    of restitution.

SourceFiles
    esBGKCollision.C

\*---------------------------------------------------------------------------*/

#ifndef esBGKCollision_H
#define esBGKCollision_H

#include "collisionKernel.H"
#include "mappedPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace collisionKernels
{

/*---------------------------------------------------------------------------*\
                    Class esBGKCollision Declaration
\*---------------------------------------------------------------------------*/

class esBGKCollision
:
    public collisionKernel
{
    // Private data

        //- Coefficient of restitution
        scalar e_;

        //- ES-BGK parameter
        scalar b_;

        //- Coefficients
        scalar a1_, b1_;

        //- Granular temperature
        volScalarField Theta_;

        //- Reference to radial distribution
//         autoPtr<radialDistributionModel> radialModel_;

        //- Coefficient related to Prandtl number
        scalar zeta_;

        //- Diameter
        tmp<volScalarField> dp_;


protected:
    // Protected data

        //- Equilibrium moments (fields)
        mappedPtrList<volScalarField> Meqf_;

        //- Equilibrium moments (fields)
        mappedList<scalar> Meq_;


public:

        //- Runtime type information
        TypeName("esBGK");


    // Constructors

        //- Construct from components
        esBGKCollision
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature,
            const bool ode
        );


    //- Destructor
    virtual ~esBGKCollision();


    // Member Functions

        //- Update equilibrium moments
        virtual void updateCells(const label celli);

        //- Return explicit collision source term
        virtual scalar explicitCollisionSource
        (
            const labelList& momentOrder,
            const label celli
        ) const;

        //- Return implicit collision source matrix
        virtual tmp<fvScalarMatrix> implicitCollisionSource
        (
            const volVelocityMoment& m
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace collisionKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
