/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 Jeffrey Heylmun
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
    Foam::ButcherTables::Euler

Description
    1st order ode time integrater

SourceFiles
    Euler.C

\*---------------------------------------------------------------------------*/

#ifndef Euler_H
#define Euler_H

#include "ButcherTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ButcherTables
{

/*---------------------------------------------------------------------------*\
                        Class Euler Declaration
\*---------------------------------------------------------------------------*/

class Euler : public ButcherTable
{
public:

    //- Runtime type information
    TypeName("Euler");

    // Constructor
    Euler(const fvMesh& mesh);

    //- Destructor
    virtual ~Euler();


    // Public Functions

        //- Number of steps
        virtual label nSteps() const
        {
            return 1;
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
