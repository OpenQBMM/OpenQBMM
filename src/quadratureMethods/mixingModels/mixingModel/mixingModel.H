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
    Foam::mixingModel

Description
    Abstract class for the run-time selection of mixing models.

SourceFiles
    mixingModel.C

\*---------------------------------------------------------------------------*/

#ifndef mixingModel_H
#define mixingModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class mixingModel Declaration
\*---------------------------------------------------------------------------*/

class mixingModel
{
    // Private data

        //- Name of the mixingModel
        const word name_;

protected:

    // Protected data

        //- Fluid-phase face-normal velocity
        const surfaceScalarField& phi_;


public:

    //- Runtime type information
    TypeName("mixingModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        mixingModel,
        dictionary,
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        ),
        (name, dict, phi)
    );


    // Constructors

        //- Construct from components
        mixingModel
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construct
        mixingModel(const mixingModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<mixingModel> New
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~mixingModel();


    // Member Functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() = 0;

        //- Solve mixing model
        virtual void solve() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mixingModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
