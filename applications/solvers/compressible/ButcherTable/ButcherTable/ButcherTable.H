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
    Foam::ButcherTable

Description
    Basic class to return Butcher tables for ode solvers.

SourceFiles
    ButcherTable.C

\*---------------------------------------------------------------------------*/

#ifndef ButcherTable_H
#define ButcherTable_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fluxIntegrator Declaration
\*---------------------------------------------------------------------------*/

class ButcherTable
{

public:

    //- Runtime type information
    TypeName("ButcherTable");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            ButcherTable,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructor
    ButcherTable(const fvMesh& mesh);

    // Selector
    static autoPtr<ButcherTable> New(const fvMesh& mesh);

    //- Destructor
    virtual ~ButcherTable();


    // Member Functions

        //- Number of steps
        virtual label nSteps() const = 0;

        //- Conserved variable coeffients
        virtual List<scalarList> conservedVariablesCoeffs() const = 0;

        //- Flux coefficients
        virtual List<scalarList> fluxCoeffs() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
