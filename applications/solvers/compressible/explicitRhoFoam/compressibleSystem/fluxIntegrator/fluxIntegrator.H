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
    Foam::fluxIntegrator

Description
    Class used to integrate the conservative fluxes using an arbirary set of
    ode coefficients defined by a Butcher table

SourceFiles
    fluxIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef fluxIntegrator_H
#define fluxIntegrator_H

#include "IOdictionary.H"
#include "ButcherTable.H"
#include "compressibleSystem.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fluxIntegrator Declaration
\*---------------------------------------------------------------------------*/

class fluxIntegrator
{
protected:

    // Protected data

        //- Name of compressible phase ('\0' by default)
        compressibleSystem& fluid_;

        //- List of ode coefficients
        autoPtr<ButcherTable> butcherTable_;

        //- Return fields that need to be set
        void setCoeffs(boolList& setFields, boolList& setDeltas);

public:

    // Constructor
    fluxIntegrator(compressibleSystem& fluid);

    //- Destructor
    virtual ~fluxIntegrator();


    // Member Functions

        //- Transport conservative variables
        virtual void integrateFluxes(const dimensionedVector& g);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
