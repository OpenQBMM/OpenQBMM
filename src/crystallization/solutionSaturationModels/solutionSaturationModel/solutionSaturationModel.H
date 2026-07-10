/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::solutionSaturationModel

Description

SourceFiles
    solutionSaturationModel.C

\*---------------------------------------------------------------------------*/

#ifndef solutionSaturationModel_H
#define solutionSaturationModel_H

#include "IOdictionary.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class solutionSaturationModel Declaration
\*---------------------------------------------------------------------------*/

class solutionSaturationModel
:
    public IOdictionary
{
    // Private Member Functions

        //- No copy construct
        solutionSaturationModel(const solutionSaturationModel&) = delete;

        //- No copy assignment
        void operator=(const solutionSaturationModel&) = delete;


public:

    //- Runtime type information
    TypeName("solutionSaturationModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        solutionSaturationModel,
        dictionary,
        (
            const dictionary& dict, const objectRegistry& db
        ),
        (dict, db)
    );


    // Constructors

        //- Default constructed (from registry)
        explicit solutionSaturationModel(const objectRegistry& db);


    // Selectors

        //- Select default constructed
        static autoPtr<solutionSaturationModel> New
        (
            const dictionary& dict,
            const objectRegistry& db
        );


    //- Destructor
    virtual ~solutionSaturationModel() = default;


    // Member Functions

        //- Saturation concentration  w.r.t. temperature
        virtual tmp<volScalarField> Csat
        (
            const volScalarField& T
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
