/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Junjie Li
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
    Foam::populationBalanceSubModels::aggregationKernels::crystalAggregation

Description
    Combined aggregation kernel for cooling crystallization that uses
    separate collision frequency and sticking efficiency sub-models.

    The aggregation kernel is calculated as:
    \f[
        K_a = C_a \cdot \beta(L_i, L_j) \cdot \psi(L_i, L_j)
    \f]

    where \f$\beta\f$ is the collision frequency and \f$\psi\f$ is the
    sticking efficiency (probability of aggregation upon collision).

    Available collision frequency models:
    - turbulent: Turbulence-induced collision
    - shear: Shear-induced collision

    Available efficiency models:
    - Hounslow: Based on crystal growth rate and shear
    - Ilievski: Based on impeller tip speed and crystal size

SourceFiles
    crystalAggregation.C

\*---------------------------------------------------------------------------*/

#ifndef crystalAggregation_H
#define crystalAggregation_H

#include "aggregationKernel.H"
#include "crystalCollisionFrequency.H"
#include "crystalAggregationEfficiency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class crystalAggregation Declaration
\*---------------------------------------------------------------------------*/

class crystalAggregation
:
    public aggregationKernel
{
    // Private data

        //- Collision frequency model
        autoPtr<crystalCollisionFrequency> frequency_;

        //- Aggregation efficiency model
        autoPtr<crystalAggregationEfficiency> efficiency_;


public:

    //- Runtime type information
    TypeName("crystalAggregation");


    // Constructors

        //- Construct from components
        crystalAggregation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~crystalAggregation();


    // Member Functions

        //- Return true if source is only dependent on size
        virtual bool pureSize() const
        {
            return false;
        }

        //- Update fields before ode solve
        virtual void preUpdate();

        //- Aggregation kernel
        virtual scalar Ka
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli,
            const label environment = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
