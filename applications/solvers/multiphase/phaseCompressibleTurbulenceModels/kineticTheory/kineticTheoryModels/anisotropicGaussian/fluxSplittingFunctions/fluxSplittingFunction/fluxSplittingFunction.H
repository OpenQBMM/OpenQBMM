/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
     \\/     M anipulation  |
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
    Foam::fluxSplittingFunction

Description

SourceFiles
    fluxSplittingFunction.C

\*---------------------------------------------------------------------------*/

#ifndef fluxSplittingFunction_H
#define fluxSplittingFunction_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                     Class fluxSplittingFunction Declaration
\*---------------------------------------------------------------------------*/

class fluxSplittingFunction
{
protected:

    // Protected data

        //- Dictionary used to construct
        const dictionary& dict_;

        //- Power used to determine how important the kinetic flux is
        const scalar h2Pow_;


public:

    //- Runtime type information
    TypeName("fluxSplittingFunction");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluxSplittingFunction,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        fluxSplittingFunction(const dictionary& dict);

        //- Disallow default bitwise copy construct
        fluxSplittingFunction(const fluxSplittingFunction&) = delete;


    // Selectors

        static autoPtr<fluxSplittingFunction> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~fluxSplittingFunction();


    // Member Functions

        virtual tmp<volScalarField> h2
        (
            const volScalarField& alpha1,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& rho1,
            const volScalarField& da,
            const volScalarField& ppfr,
            const dimensionedScalar& e
        ) const = 0;

        virtual bool read()
        {
            return true;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fluxSplittingFunction&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
