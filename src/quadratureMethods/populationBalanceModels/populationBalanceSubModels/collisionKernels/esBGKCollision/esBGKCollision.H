/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2018 by Alberto Passalacqua
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
    Foam::populationBalanceSubModels::collisionKernels::esBGKCollision

Description
    Extension of the BGK collisional model that includes the coefficient
    of restitution.

SourceFiles
    esBGKCollision.C

\*---------------------------------------------------------------------------*/

#ifndef esBGKCollision_H
#define esBGKCollision_H

#include "BGKCollision.H"
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
    public BGKCollision
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


protected:

    // Protected Functions

        //- Return the covariance tensor using quadrature moments
        virtual symmTensor covariance
        (
            const label celli,
            const scalar& u,
            const scalar& v,
            const scalar& w
        );

        //- Return the covariance tensor given a list of moments
        virtual symmTensor covariance
        (
            const mappedScalarList& moments,
            const scalar& u,
            const scalar& v,
            const scalar& w
        );

        //- Return the coefficient of restitution
        virtual scalar e() const
        {
            return e_;
        }


public:

        //- Runtime type information
        TypeName("esBGK");


    // Constructors

        //- Construct from components
        esBGKCollision
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const velocityQuadratureApproximation& quadrature
        );


    //- Destructor
    virtual ~esBGKCollision();


    // Member Functions

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
