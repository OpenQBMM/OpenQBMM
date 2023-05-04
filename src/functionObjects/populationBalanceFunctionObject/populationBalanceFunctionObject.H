/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022-2023 Alberto Passalacqua
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
    Foam::functionObjects::populationBalance

Description
    Solves a population balance equation with quadrature-based moment methods.
    It should be executed in association to the fluid module solver.

    Example:
    \verbatim
        pbe
        {
            libs ("libOpenQBMMFunctionObjects.so");

            type populationBalance;

            phi phi;  // optional
            rho rho;  // optional
        }
    \endverbatim

    See /test/populationBalanceModel/functionObjectTestCase for an example.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    populationBalance.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_populationBalance_H
#define functionObjects_populationBalance_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class populationBalance Declaration
\*---------------------------------------------------------------------------*/

class populationBalance
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name of field to process
        word fieldName_;

        //- Name of flux field (optional)
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;

        //- Population balance properties dictionary
        IOdictionary populationBalanceProperties_;

        //- Volumetric flux field
        surfaceScalarField phi_;

        //- Population balance model to be solved
        autoPtr<populationBalanceModel> populationBalanceModelPtr_;

    // Private Member Functions

        void updateVolumetricFaceFlux();


public:

    //- Runtime type information
    TypeName("populationBalance");


    // Constructors

        //- Construct from Time and dictionary
        populationBalance
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        populationBalance(const populationBalance&) = delete;


    //- Destructor
    virtual ~populationBalance();


    // Member Functions

        //- Read the populationBalance data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Calculate the populationBalance
        virtual bool execute();

        //- Write
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const populationBalance&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //