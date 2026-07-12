/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2023 by Your Name
    Copyright (C) 2023 Your Organization
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

SourceFiles
    aggregationTimeModel.C

\*---------------------------------------------------------------------------*/
#ifndef aggregationTimeModel_H
#define aggregationTimeModel_H

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class aggregationTimeModel Declaration
\*---------------------------------------------------------------------------*/

class aggregationTimeModel
:
    public aggregationKernel
{
    // Private data

        //- Field to store the aggregation rate
        volScalarField FField_;

        //- Reference to the local growth-rate field
        const volScalarField& growthRate_;

        //- Reference to the local supersaturation field
        const volScalarField& sigma_;

        //- Reference to the local turbulent dissipation rate
        const volScalarField& epsilon_;

        //- Reference to the local kinematic viscosity field
        const volScalarField& nu_;

        //- Reference to the local temperature field
        const volScalarField& T_;

        //- Reference liquid density used to obtain dynamic viscosity
        scalar rhoLiquid_;


public:

    //- Runtime type information
    TypeName("aggregationTimeModel");


    // Constructors

        //- Construct from components
        aggregationTimeModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~aggregationTimeModel();


    // Member Functions

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
