/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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
    Foam::bubblePressureModels::BeisheuvelBubblePressure

Description
    Wall lubrication model of Tomiyama.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
         in a Gas-Liquid System within an Euler-Euler Framework"
        Otromke, M.
        PhD Thesis, April 2013
    \endverbatim

SourceFiles
    BeisheuvelBubblePressure.C

\*---------------------------------------------------------------------------*/

#ifndef BeisheuvelBubblePressure_H
#define BeisheuvelBubblePressure_H

#include "bubblePressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace bubblePressureModels
{

/*---------------------------------------------------------------------------*\
                   Class BeisheuvelBubblePressure Declaration
\*---------------------------------------------------------------------------*/

class Beisheuvel
:
    public bubblePressureModel
{
    // Private data

        //- Characteristic channel dimension
        const dimensionedScalar Cbp1_;


public:

    //- Runtime type information
    TypeName("Beisheuvel");


    // Constructors

        //- Construct from components
        Beisheuvel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Beisheuvel();


    // Member Functions

        //- Return bubble pressure
        tmp<volScalarField> pb(const label, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace bubblePressureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
