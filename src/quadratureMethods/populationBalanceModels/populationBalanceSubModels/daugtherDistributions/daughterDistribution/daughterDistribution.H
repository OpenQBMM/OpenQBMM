/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::daughterDistribution

Description
    Abstract class for the daughter distribution of fragments in breakup
    processes.
    
SourceFiles
    daughterDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef daughterDistribution_H
#define daughterDistribution_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{

/*---------------------------------------------------------------------------*\
                    Class daughterDistribution Declaration
\*---------------------------------------------------------------------------*/

class daughterDistribution
{
    // Private member functions

        //- Disallow default bitwise copy construct
        daughterDistribution(const daughterDistribution&);

        //- Disallow default bitwise assignment
        void operator=(const daughterDistribution&);


protected:

    // Protected data

        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("daughterDistribution");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        daughterDistribution,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        daughterDistribution(const dictionary& dict);


    // Selectors

        static autoPtr<daughterDistribution> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~daughterDistribution();


    // Member Functions

        //- Moment of the daughter distribution function
        virtual tmp<volScalarField> mD
        (
            const label order,
            const volScalarField& abscissa
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
