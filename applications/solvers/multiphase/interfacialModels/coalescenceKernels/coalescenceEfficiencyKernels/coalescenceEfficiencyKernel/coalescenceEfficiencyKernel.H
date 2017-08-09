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
    Foam::populationBalanceSubModels::coalescenceEfficiencyKernel

Description
    Abstract class for coalescence effeciency kernels.

SourceFiles
    coalescenceEfficiencyKernel.C

\*---------------------------------------------------------------------------*/

#ifndef coalescenceEfficiencyKernel_H
#define coalescenceEfficiencyKernel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coalescenceEfficiencyKernel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceEfficiencyKernel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        coalescenceEfficiencyKernel(const coalescenceEfficiencyKernel&);

        //- Disallow default bitwise assignment
        void operator=(const coalescenceEfficiencyKernel&);


protected:

    // Protected data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("coalescenceEfficiencyKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coalescenceEfficiencyKernel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        coalescenceEfficiencyKernel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Selectors

        static autoPtr<coalescenceEfficiencyKernel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~coalescenceEfficiencyKernel();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<volScalarField> Pc
        (
            const label nodei,
            const label nodej
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
