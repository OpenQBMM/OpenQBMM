/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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
    Foam::mixingSubModels::diffusionModel

Description
    Abstract class for moment diffusion models.

SourceFiles
    diffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionModel_H
#define diffusionModel_H

#include "dictionary.H"
#include "volFields.H"
#include "fvMatrix.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{

/*---------------------------------------------------------------------------*\
                    Class diffusionModel Declaration
\*---------------------------------------------------------------------------*/

class diffusionModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        diffusionModel(const diffusionModel&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionModel&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("diffusionModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        diffusionModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        diffusionModel(const dictionary& dict);


    // Selectors

        static autoPtr<diffusionModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~diffusionModel();


    // Member Functions

        //- Return the source term for diffusion term in the moment equation
        virtual tmp<fvScalarMatrix> momentDiff
        (
            const volScalarField& moment
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
