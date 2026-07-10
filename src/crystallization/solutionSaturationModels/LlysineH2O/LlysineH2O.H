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
    Foam::solutionSaturationModels::LlysineH2O

Description
    LlysineH2O equation for the saturation vapour temperature in terms of
    the vapour pressure (in Pa).

    \f[
        T_sat = \sum_i C_i p^i
    \f]

    where \f$p\f$ is the pressure in Pa and \f$C\f$ are the coefficients.

    Currently this class only provides \f$T_sat\f$, the inverse function to
    return the vapour pressure for a given temperature are not implemented.

SourceFiles
    LlysineH2O.C

\*---------------------------------------------------------------------------*/

#ifndef solutionSaturationModels_LlysineH2O_H
#define solutionSaturationModels_LlysineH2O_H

#include "solutionSaturationModel.H"
#include "LlysineH2O.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solutionSaturationModels
{

/*---------------------------------------------------------------------------*\
                         Class LlysineH2O Declaration
\*---------------------------------------------------------------------------*/

class LlysineH2O
:
    public solutionSaturationModel
{
    // Private data


public:

    //- Runtime type information
    TypeName("LlysineH2O");

    // Constructors

        //- Construct from a dictionary
        LlysineH2O(const dictionary& dict, const objectRegistry& db);


    //- Destructor
    virtual ~LlysineH2O();


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
