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
    Foam::solutionSaturationModels::KclH2o

Description
    KclH2o equation for the saturation vapour temperature in terms of
    the vapour pressure (in Pa).

    \f[
        T_sat = \sum_i C_i p^i
    \f]

    where \f$p\f$ is the pressure in Pa and \f$C\f$ are the coefficients.

    Currently this class only provides \f$T_sat\f$, the inverse function to
    return the vapour pressure for a given temperature are not implemented.

SourceFiles
    KclH2o.C

\*---------------------------------------------------------------------------*/

#ifndef solutionSaturationModels_KclH2o_H
#define solutionSaturationModels_KclH2o_H

#include "solutionSaturationModel.H"
#include "KclH2o.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solutionSaturationModels
{

/*---------------------------------------------------------------------------*\
                         Class KclH2o Declaration
\*---------------------------------------------------------------------------*/

class KclH2o
:
    public solutionSaturationModel
{
    // Private data


public:

    //- Runtime type information
    TypeName("KclH2o");

    // Constructors

        //- Construct from a dictionary
        KclH2o(const dictionary& dict, const objectRegistry& db);


    //- Destructor
    virtual ~KclH2o();


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
