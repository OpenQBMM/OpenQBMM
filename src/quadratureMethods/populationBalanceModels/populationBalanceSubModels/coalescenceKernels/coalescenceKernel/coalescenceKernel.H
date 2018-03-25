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
    Foam::populationBalanceSubModels::coalescenceKernel

Description
    Abstract class for coalescence kernels.
    
SourceFiles
    coalescenceKernel.C

\*---------------------------------------------------------------------------*/

#ifndef coalescenceKernel_H
#define coalescenceKernel_H

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
                    Class coalescenceKernel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceKernel
{   
    // Private member functions

        //- Disallow default bitwise copy construct
        coalescenceKernel(const coalescenceKernel&);

        //- Disallow default bitwise assignment
        void operator=(const coalescenceKernel&);


protected:

    // Protected data

        const dictionary& dict_;

        //- Coefficient of coalescence kernel
        const dimensionedScalar Cc_;


public:

    //- Runtime type information
    TypeName("coalescenceKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coalescenceKernel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        coalescenceKernel(const dictionary& dict);


    // Selectors

        static autoPtr<coalescenceKernel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~coalescenceKernel();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<volScalarField> Ka
        (
            const volScalarField& abscissa1,
            const volScalarField& abscissa2
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
