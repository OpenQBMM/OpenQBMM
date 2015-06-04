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
    Foam::populationBalanceSubModels::breakupKernel

Description
    Abstract class for breakup kernels.
    
SourceFiles
    breakupKernel.C

\*---------------------------------------------------------------------------*/

#ifndef breakupKernel_H
#define breakupKernel_H

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
                    Class breakupKernel Declaration
\*---------------------------------------------------------------------------*/

class breakupKernel
{   
    // Private member functions

        //- Disallow default bitwise copy construct
        breakupKernel(const breakupKernel&);

        //- Disallow default bitwise assignment
        void operator=(const breakupKernel&);


protected:

    // Protected data

        const dictionary& dict_;

        //- Coefficient of aggregation kernel
        const dimensionedScalar Cb_;


public:

    //- Runtime type information
    TypeName("breakupKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        breakupKernel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        breakupKernel(const dictionary& dict);


    // Selectors

        static autoPtr<breakupKernel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~breakupKernel();


    // Member Functions

        //- Breakup kernel
        virtual tmp<volScalarField> Kb
        (
            const volScalarField& abscissa
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
