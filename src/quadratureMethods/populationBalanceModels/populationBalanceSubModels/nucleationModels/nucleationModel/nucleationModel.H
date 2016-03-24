/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alberto Passalacqua
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
    // Private member functions

        //- Disallow default bitwise copy construct
        nucleationModel(const nucleationModel&);

        //- Disallow default bitwise assignment
        void operator=(const nucleationModel&);


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
        virtual tmp<volScalarField> nucleationSource
        (
            const volUnivariateMoment& moment
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
