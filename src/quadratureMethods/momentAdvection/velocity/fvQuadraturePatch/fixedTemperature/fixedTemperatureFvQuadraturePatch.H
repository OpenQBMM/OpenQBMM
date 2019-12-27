/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 Alberto Passalacqua
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
    Foam::fixedTemperatureFvQuadraturePatch

Description
    Abstract class for the advection of a velocity moment set.

SourceFiles
    fixedTemperatureFvQuadraturePatch.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTemperatureFvQuadraturePatch_H
#define fixedTemperatureFvQuadraturePatch_H

#include "fvQuadraturePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class fixedTemperatureFvQuadraturePatch Declaration
\*---------------------------------------------------------------------------*/

class fixedTemperatureFvQuadraturePatch
:
    public fvQuadraturePatch
{
    // Private data

        //- Fixed wall temperature
        scalarField wallTemperature_;

        //- Number of veloicity components
        const label nVelocityCmpts_;

        //- Moment orders
        labelList order000_;
        labelList order100_;
        labelList order010_;
        labelList order001_;
        labelList order200_;
        labelList order020_;
        labelList order002_;

public:

    //- Runtime type information
    TypeName("fixedTemperature");


    // Constructors

        //- Construct from velocityMomentSet
        fixedTemperatureFvQuadraturePatch
        (
            const fvPatch& patch,
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            PtrList<surfaceVelocityNode>& nodesOwn,
            PtrList<surfaceVelocityNode>& nodesNei
        );


    //- Destructor
    virtual ~fixedTemperatureFvQuadraturePatch();

    // Member Functions

        //- Update boundary flux
        virtual void update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
