/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::coalesenceEffeciencyKernel

Description
    Abstract class for coalesence effeciency kernels.

SourceFiles
    coalesenceEffeciencyKernel.C

\*---------------------------------------------------------------------------*/

#ifndef coalesenceEffeciencyKernel_H
#define coalesenceEffeciencyKernel_H

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
                    Class coalesenceEffeciencyKernel Declaration
\*---------------------------------------------------------------------------*/

class coalesenceEffeciencyKernel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        coalesenceEffeciencyKernel(const coalesenceEffeciencyKernel&);

        //- Disallow default bitwise assignment
        void operator=(const coalesenceEffeciencyKernel&);


protected:

    // Protected data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("coalesenceEffeciencyKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coalesenceEffeciencyKernel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        coalesenceEffeciencyKernel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Selectors

        static autoPtr<coalesenceEffeciencyKernel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~coalesenceEffeciencyKernel();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<volScalarField> Pc
        (
            const volScalarField& d1,
            const volScalarField& d2
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
