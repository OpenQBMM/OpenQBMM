/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::bubbleBreakupKernel

Description
    Class to calculate the breakup rate for bubbles

SourceFiles
    bubbleBreakupKernel.C

\*---------------------------------------------------------------------------*/

#ifndef bubbleBreakupKernel_H
#define bubbleBreakupKernel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class bubbleBreakupKernel Declaration
\*---------------------------------------------------------------------------*/

class bubbleBreakupKernel
{
protected:
    // Protected member functions

        //- Disallow default bitwise copy construct
        bubbleBreakupKernel(const bubbleBreakupKernel&);

        //- Disallow default bitwise assignment
        void operator=(const bubbleBreakupKernel&);


    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Coefficient of beakup kernel
        const dimensionedScalar Cb_;


public:

    //- Runtime type information
    TypeName("bubbleBreakupKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        bubbleBreakupKernel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        bubbleBreakupKernel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Selectors

        static autoPtr<bubbleBreakupKernel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~bubbleBreakupKernel();


    // Member Functions

        //- Breakup kernel
        virtual tmp<volScalarField> Kb(const label nodei) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
