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

        //- Diameter
        dimensionedScalar dp_;

        //- Coefficient related to Prandtl number
        scalar zeta_;


protected:
    // Protected data

        //- Equilibrium moments (fields)
        mappedPtrList<volScalarField> Meqf_;

        //- Equilibrium moments (fields)
        mappedList<scalar> Meq_;


    // Protected functions

        //- Update equilibrium moments for a cell (1D)
        void updateCells1D(const label celli);

        //- Update equilibrium moments for a cell (2D)
        void updateCells2D(const label celli);

        //- Update equilibrium moments for a cell (3D)
        void updateCells3D(const label celli);

        //- Update equilibrium moments for fields (1D)
        void updateFields1D();

        //- Update equilibrium moments for fields (2D)
        void updateFields2D();

        //- Update equilibrium moments for fields (3D)
        void updateFields3D();


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

        //- Update equilibrium moments
        virtual void updateFields();

        //- Return explicit collision source term
        virtual scalar explicitCollisionSource
        (
            const label mi,
            const label celli
        ) const;

        //- Return implicit collision source matrix
        virtual tmp<fvScalarMatrix> implicitCollisionSource
        (
            const volVectorMoment& m
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace collisionKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
