/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2023 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::nucleationModel

Description
    Abstract class for nucleation model.

SourceFiles
    nucleationModel.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationModel_H
#define nucleationModel_H

#include "dictionary.H"
#include "volFields.H"
#include "moments.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{

/*---------------------------------------------------------------------------*\
                    Class nucleationModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationModel
{
protected:

    // Protected data

        const dictionary& dict_;

        //- Reference to mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("nucleationModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        nucleationModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        nucleationModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        nucleationModel(const nucleationModel&) = delete;


    // Selectors

        static autoPtr<nucleationModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~nucleationModel();


    // Member Functions

        //- Nucleation rate
        virtual scalar nucleationSource
        (
            const label& momentOrder,
            const label celli,
            const label environment = 0
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const nucleationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
