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
    Foam::populationBalanceModel

Description
    Abstract class for the run-time selection of the population balance model.

SourceFiles
    populationBalanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class populationBalanceModel Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceModel
{
    // Private data

        //- Name of the populationBalanceModel
        const word name_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        populationBalanceModel(const populationBalanceModel&);

        //- Disallow default bitwise assignment
        void operator=(const populationBalanceModel&);


protected:

    // Protected data

        //- Fluid-phase face-normal velocity
        const surfaceScalarField& phi_;


public:

    //- Runtime type information
    TypeName("populationBalanceModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        populationBalanceModel,
        dictionary,
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        ),
        (name, dict, U, phi)
    );


    // Constructors

        //- Construct from components
        populationBalanceModel
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Selectors

        //- Select null constructed
        static autoPtr<populationBalanceModel> New
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~populationBalanceModel();


    // Member Functions

        //- Solve population balance equation
        virtual void solve() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
