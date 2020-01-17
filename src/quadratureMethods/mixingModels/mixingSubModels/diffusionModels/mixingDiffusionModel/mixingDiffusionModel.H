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
    Foam::mixingSubModels::mixingDiffusionModel

Description
    Abstract class for moment diffusion models.

SourceFiles
    diffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef mixingDiffusionModel_H
#define mixingDiffusionModel_H

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
                    Class mixingDiffusionModel Declaration
\*---------------------------------------------------------------------------*/

class mixingDiffusionModel
{
protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("mixingDiffusionModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        mixingDiffusionModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        mixingDiffusionModel(const dictionary& dict);

        //- Disallow default bitwise copy construct
        mixingDiffusionModel(const mixingDiffusionModel&) = delete;


    // Selectors

        static autoPtr<mixingDiffusionModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~mixingDiffusionModel();


    // Member Functions

        //- Return the source term for diffusion term in the moment equation
        virtual tmp<fvScalarMatrix> momentDiff
        (
            const volScalarField& moment
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mixingDiffusionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
