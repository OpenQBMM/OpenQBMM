/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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
    Foam::mixingSubModels::mixingKernels::FokkerPlanck

Description
    Fokker-Planck mixing model.

    References
    \verbatim
        "Computational Models for Turbulent Reacting Flows"
        R O Fox
        Cambridge University Press
    \endverbatim

SourceFiles
    FokkerPlanck.C

\*---------------------------------------------------------------------------*/

#ifndef FokkerPlanck_H
#define FokkerPlanck_H

#include "mixingKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace mixingKernels
{

/*---------------------------------------------------------------------------*\
                        Class FokkerPlanck Declaration
\*---------------------------------------------------------------------------*/

class FokkerPlanck
:
    public mixingKernel
{
public:

        //- Runtime type information
        TypeName("FokkerPlanck");


    // Constructors

        //- Construct from components
        FokkerPlanck
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~FokkerPlanck();


    // Member Functions

        //- Mixing kernel
        virtual tmp<fvScalarMatrix> K
        (
            const volScalarMoment& moment,
            const volScalarMomentFieldSet& moments
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingKernels
} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
