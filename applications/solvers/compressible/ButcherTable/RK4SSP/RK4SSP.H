/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 Jeffrey Heylmun
    Copyright (C) 2021 Alberto Passalacqua
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
        Foam::ButcherTables::RK3

Description
    Second order Runge-Kutta ode time integrater

SourceFiles
    RK4.C

\*---------------------------------------------------------------------------*/

#ifndef RK4SSP_H
#define RK4SSP_H

#include "ButcherTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ButcherTables
{

/*---------------------------------------------------------------------------*\
                        Class RK4 Declaration
\*---------------------------------------------------------------------------*/

class RK4SSP : public ButcherTable
{
public:

    //- Runtime type information
    TypeName("RK4SSP");

    // Constructor
    RK4SSP(const fvMesh& mesh);

    //- Destructor
    virtual ~RK4SSP();

    // Public Functions

        //- Number of steps
        virtual label nSteps() const
        {
            return 10;
        }

        //- Conserved variable coeffients
        virtual List<scalarList> conservedVariablesCoeffs() const;

        //- Flux coefficients
        virtual List<scalarList> fluxCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ButcherTables
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
