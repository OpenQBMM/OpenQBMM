/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    Foam::momentGenerationSubModels::alphaAndDiameterVelocity

Description
    Moment generation done using volume fraction, density, bubble diameter,
    and velocity.

SourceFiles
    alphaAndDiameterVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef alphaAndDiameterVelocity_H
#define alphaAndDiameterVelocity_H

#include "alphaAndDiameter.H"
#include "rhoThermo.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{

/*---------------------------------------------------------------------------*\
                    Class alphaAndDiameterVelocity Declaration
\*---------------------------------------------------------------------------*/

class alphaAndDiameterVelocity
:
    public alphaAndDiameter
{
public:

        //- Runtime type information
        TypeName("alphaAndDiameterVelocity");


    // Constructors

        //- Construct from components
        alphaAndDiameterVelocity
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const labelListList& momentsOrders,
            const label nNodes
        );


    //- Destructor
    virtual ~alphaAndDiameterVelocity();


    // Public member functions

        //- Update weights and abscissae given a dictionary
        virtual void updateMoments
        (
            const dictionary& dict,
            const label patchi = -1
        );

        //- Update weights and abscissae given a dictionary and cell set
        virtual void updateMoments
        (
            const dictionary& dict,
            const labelList& cells
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace momentGenerationSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
