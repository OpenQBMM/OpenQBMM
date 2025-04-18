/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2023 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::coalescenceFrequencyKernel

Description
    Abstract class for coalescence effeciency kernels.

SourceFiles
    coalescenceFrequencyKernel.C

\*---------------------------------------------------------------------------*/

#ifndef coalescenceFrequencyKernel_H
#define coalescenceFrequencyKernel_H

#include "turbulentFluidThermoModel.H"
#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class coalescenceFrequencyKernel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceFrequencyKernel
{
public:

    //- Runtime type information
    TypeName("coalescenceFrequencyKernel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        coalescenceFrequencyKernel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& continuousPhase
        ),
        (dict, mesh, continuousPhase)
    );


    // Constructors

        //- Construct from components
        coalescenceFrequencyKernel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& continuousPhase
        );

        //- Disallow default bitwise copy construct
        coalescenceFrequencyKernel(const coalescenceFrequencyKernel&) = delete;


    // Selectors

        static autoPtr<coalescenceFrequencyKernel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& continuousPhase
        );


    //- Destructor
    virtual ~coalescenceFrequencyKernel();


    // Member Functions

        //- Update fields
        virtual void update
        (
            const fluidThermo& thermo,
            const turbulenceModel& turb
        )
        {}

        //- Aggregation kernel
        virtual scalar omega
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const coalescenceFrequencyKernel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
