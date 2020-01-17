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
    Foam::populationBalanceSubModels::aggregationKernel

Description
    Abstract class for aggregation kernels.

SourceFiles
    aggregationKernel.C

\*---------------------------------------------------------------------------*/

#ifndef aggregationKernel_H
#define aggregationKernel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{

/*---------------------------------------------------------------------------*\
                    Class aggregationKernel Declaration
\*---------------------------------------------------------------------------*/

class aggregationKernel
{
protected:

    // Protected data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Coefficient of aggregation kernel
        const dimensionedScalar Ca_;


    // Protected functions

        //- Return true if source is only depenent on size
        virtual bool pureSize() const
        {
            return true;
        }

        //- Length based aggregation source given two abscissae and an order
        scalar nodeSource
        (
            const scalar& abscissa1,
            const scalar& abscissa2,
            const label momentOrder
        ) const;

        //- Mass based aggregation source given two abscissae and an order
        scalar massNodeSource
        (
            const scalar& abscissa1,
            const scalar& abscissa2,
            const label momentOrder
        ) const;


public:

    //- Runtime type information
    TypeName("aggregationKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        aggregationKernel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        aggregationKernel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        aggregationKernel(const aggregationKernel&) = delete;


    // Selectors

        static autoPtr<aggregationKernel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~aggregationKernel();


    // Member Functions

        //- Update fields before ode solve
        virtual void preUpdate()
        {
            return;
        }

        //- Aggregation kernel
        virtual scalar Ka
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli,
            const label environment = 0
        ) const = 0;

        //- Return aggregation source terms
        scalar aggregationSource
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature,
            const label enviroment
        );

        //- Return aggregation source terms
        scalar aggregationSource
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature,
            const label enviroment
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const aggregationKernel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
