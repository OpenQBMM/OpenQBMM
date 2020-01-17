/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2017-2018 by Alberto Passalacqua
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
    Foam::populationBalanceSubModels::aggregationKernels::coalescence

Description
    Class to calculate the coalescence rate for bubbles based on coalescence
    frequency and efficiency.

SourceFiles
    coalescenceKernel.C

\*---------------------------------------------------------------------------*/

#ifndef coalescenceKernel_H
#define coalescenceKernel_H

#include "aggregationKernel.H"
#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "turbulentFluidThermoModel.H"
#include "coalescenceFrequencyKernel.H"
#include "coalescenceEfficiencyKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class coalescenceKernel Declaration
\*---------------------------------------------------------------------------*/

class coalescence
:
    public aggregationKernel
{
    // Private data

        //- continuous phase name
        const word continuousPhase_;

        //- Coalescence frequency model
        autoPtr<coalescenceFrequencyKernel> frequency_;

        //- Coalescence effeciency model
        autoPtr<coalescenceEfficiencyKernel> efficiency_;


public:

        //- Runtime type information
        TypeName("coalescence");

    // Constructor
    coalescence
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~coalescence();

    // Member Functions

        //- Return true if source is only depenent on size
        virtual bool pureSize() const
        {
            return false;
        }

        //- Update fields
        virtual void preUpdate();

        //- Coalescence kernel for a cell
        scalar Ka
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli,
            const label environment = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
