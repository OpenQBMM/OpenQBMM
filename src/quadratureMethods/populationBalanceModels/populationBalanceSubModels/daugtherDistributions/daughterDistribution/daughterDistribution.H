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

        //- Disallow default bitwise copy construct
        daughterDistribution(const daughterDistribution&) = delete;


    // Selectors

        static autoPtr<daughterDistribution> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~daughterDistribution();


    // Member Functions

        //- Moment of the daughter distribution function (length based)
        virtual scalar mD
        (
            const label& order,
            const scalar& abscissa
        ) const = 0;

        //- Moment of the daughter distribution function (mass/volume based)
        virtual scalar mDMass
        (
            const label& order,
            const scalar& abscissa
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const daughterDistribution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
