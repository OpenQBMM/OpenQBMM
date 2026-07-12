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
    Foam::populationBalanceSubModels::aggregationKernels::
    crystalCollisionFrequencies::shearFrequency

Description
    Shear-induced collision frequency model for crystal aggregation.

    The collision frequency is calculated as:
    \f[
        \beta_{shear}(L_i, L_j) = \frac{1}{6} \dot{\gamma} (L_i + L_j)^3
    \f]

    \vartable
        \dot{\gamma} | Shear rate [1/s]
        L_i, L_j     | Crystal sizes [m]
    \endvartable

SourceFiles
    shearFrequency.C

\*---------------------------------------------------------------------------*/

#ifndef shearFrequency_H
#define shearFrequency_H

#include "crystalCollisionFrequency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalCollisionFrequencies
{

/*---------------------------------------------------------------------------*\
                      Class shearFrequency Declaration
\*---------------------------------------------------------------------------*/

class shearFrequency
:
    public crystalCollisionFrequency
{
    // Private data

        //- Shear coefficient (default 1/6)
        scalar Cshear_;

        //- Reference to kinematic viscosity field
        const volScalarField& nu_;

        //- Reference to turbulent dissipation rate
        const volScalarField& epsilon_;


public:

    //- Runtime type information
    TypeName("shear");


    // Constructors

        //- Construct from components
        shearFrequency
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~shearFrequency();


    // Member Functions

        //- Calculate collision frequency
        virtual scalar beta
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace crystalCollisionFrequencies
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
