/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-06-26 Jeff Heylmun:    Changed alpha to phase so that twoPhaseSystem can
                            be accessed
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::kineticTheoryModels::fluxSplittingFunctions::alphG0

Description

SourceFiles
    particlePressure.C

\*---------------------------------------------------------------------------*/

#ifndef particlePressure_H
#define particlePressure_H

#include "fluxSplittingFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace fluxSplittingFunctions
{

/*---------------------------------------------------------------------------*\
                        Class particlePressure Declaration
\*---------------------------------------------------------------------------*/

class particlePressure
:
    public fluxSplittingFunction
{
    // Private data

        //- Minimum pressure
        const dimensionedScalar minPpk_;


public:

    //- Runtime type information
    TypeName("particlePressure");


    // Constructors

        //- Construct from components
        particlePressure(const dictionary& dict);


    //- Destructor
    virtual ~particlePressure();


    // Member functions

        tmp<volScalarField> h2
        (
            const volScalarField& alpha1,
            const volScalarField& Theta,
            const volScalarField& g0,
            const volScalarField& rho1,
            const volScalarField& da,
            const volScalarField& ppfr,
            const dimensionedScalar& e
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluxSplittingFunctions
} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
