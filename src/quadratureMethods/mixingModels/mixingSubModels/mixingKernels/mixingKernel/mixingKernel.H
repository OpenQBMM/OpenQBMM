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
    Foam::mixingSubModels::mixingKernel

Description
    Abstract class for mixing model kernels.

SourceFiles
    mixingKernel.C

\*---------------------------------------------------------------------------*/

#ifndef mixingKernel_H
#define mixingKernel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "momentFieldSets.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{

/*---------------------------------------------------------------------------*\
                    Class mixingKernel Declaration
\*---------------------------------------------------------------------------*/

class mixingKernel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        mixingKernel(const mixingKernel&);

        //- Disallow default bitwise assignment
        void operator=(const mixingKernel&);


protected:

    // Protected data

        const dictionary& dict_;

        //- Coefficient of aggregation kernel
        const dimensionedScalar Cphi_;

        //- Coefficient of relaxation control
        const dimensionedScalar Cmixing_;


public:

    //- Runtime type information
    TypeName("mixingKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        mixingKernel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        mixingKernel(const dictionary& dict);


    // Selectors

        static autoPtr<mixingKernel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~mixingKernel();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<fvScalarMatrix> K
        (
            const volMoment& moment,
            const volMomentFieldSet& moments
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
