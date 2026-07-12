/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::solutionSaturationModels::K2SO4H2O

Description
    K2SO4H2O temperature-dependent solubility model.

    \f[
        C^*(T) = -249.369 + 0.3761 T + 0.0029415 T^2
    \f]

    where \f$T\f$ is in K and \f$C^*\f$ is in kg/m3.

SourceFiles
    K2SO4H2O.C

\*---------------------------------------------------------------------------*/

#ifndef solutionSaturationModels_K2SO4H2O_H
#define solutionSaturationModels_K2SO4H2O_H

#include "solutionSaturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solutionSaturationModels
{

/*---------------------------------------------------------------------------*\
                         Class K2SO4H2O Declaration
\*---------------------------------------------------------------------------*/

class K2SO4H2O
:
    public solutionSaturationModel
{
    // Private data


public:

    //- Runtime type information
    TypeName("K2SO4H2O");

    // Constructors

        //- Construct from a dictionary
        K2SO4H2O(const dictionary& dict, const objectRegistry& db);


    //- Destructor
    virtual ~K2SO4H2O();


    // Member Functions

        //- Saturation concentration  w.r.t. temperature
        virtual tmp<volScalarField> Csat
        (
            const volScalarField& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solutionSaturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
