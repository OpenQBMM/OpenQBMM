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
    Foam::populationBalanceModel

Description
    Abstract class for the run-time selection of the population balance model.

SourceFiles
    populationBalanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "IOdictionary.H"
#include "fvCFD.H"
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
:
    public regIOobject
{
    // Private data

        //- Name of the populationBalanceModel
        const word name_;


protected:

    // Protected data

        //- Reference to populationBalanceProperties
        //  Non-const so ode solver can read on modification
        IOdictionary& populationBalanceProperties_;

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
            const surfaceScalarField& phi
        ),
        (name, dict, phi)
    );


    // Constructors

        //- Construct from components
        populationBalanceModel
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construct
        populationBalanceModel(const populationBalanceModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<populationBalanceModel> New
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~populationBalanceModel();


    // Member Functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const = 0;

        //- Return the courant number based on abscissae
        virtual scalar CoNum() const = 0;

        //- Solve population balance equation
        virtual void solve() = 0;


    // regIOobject functions

        //- Read
        virtual bool readIfModified()
        {
            return true;
        }

        //- Read
        virtual bool writeData(Ostream&) const
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const populationBalanceModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
