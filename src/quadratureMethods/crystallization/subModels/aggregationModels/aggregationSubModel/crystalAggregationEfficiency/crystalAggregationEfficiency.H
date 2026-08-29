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
    Foam::populationBalanceSubModels::crystalAggregationEfficiency

Description
    Abstract base class for crystal aggregation efficiency models in 
    cooling crystallization processes.

    The efficiency represents the probability that two colliding crystals
    will stick together and form an aggregate.

SourceFiles
    crystalAggregationEfficiency.C
    newCrystalAggregationEfficiency.C

\*---------------------------------------------------------------------------*/

#ifndef crystalAggregationEfficiency_H
#define crystalAggregationEfficiency_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                Class crystalAggregationEfficiency Declaration
\*---------------------------------------------------------------------------*/

class crystalAggregationEfficiency
{
protected:

    // Protected data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("crystalAggregationEfficiency");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        crystalAggregationEfficiency,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        crystalAggregationEfficiency
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        crystalAggregationEfficiency
        (
            const crystalAggregationEfficiency&
        ) = delete;


    // Selectors

        static autoPtr<crystalAggregationEfficiency> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~crystalAggregationEfficiency();


    // Member Functions

        //- Pre-update hook for field updates
        virtual void update()
        {}

        //- Calculate aggregation efficiency (sticking probability)
        //  Returns value between 0 and 1
        virtual scalar Pc
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const crystalAggregationEfficiency&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
