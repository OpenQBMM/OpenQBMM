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
    Foam::populationBalanceSubModels::coalescenceEfficiencyKernels::
    constant

Description
    Constant coalescence efficiency

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "coalescenceEfficiencyKernel.H"
#include "twoPhaseSystem.H"
#include "phaseModel.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceEfficiencyKernels
{

/*---------------------------------------------------------------------------*\
                    Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public coalescenceEfficiencyKernel
{
    // Private data

        //- Refrence to two phase system
        const twoPhaseSystem& fluid_;

        //- Efficiency coefficient
        dimensionedScalar Ceff_;


public:

        //- Runtime type information
        TypeName("constant");


    // Constructors

        //- Construct from components
        constant
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<volScalarField> Pc
        (
            const label nodei,
            const label nodej
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceEfficiencyKernels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
