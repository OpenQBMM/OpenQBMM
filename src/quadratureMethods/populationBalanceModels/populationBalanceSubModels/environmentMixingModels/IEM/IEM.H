/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2016-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::environmentMixingModel::IEM

Description
    Interaction by exchange with the mean mixing model.

    References
    \verbatim
        "Computational Models for Turbulent Reacting Flows"
        R O Fox
        Cambridge University Press
    \endverbatim

SourceFiles
    IEM.C

\*---------------------------------------------------------------------------*/

#ifndef IEM_H
#define IEM_H

#include "environmentMixingModel.H"
#include "turbulentFluidThermoModel.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace environmentMixingModels
{

/*---------------------------------------------------------------------------*\
                           Class IEM Declaration
\*---------------------------------------------------------------------------*/

class IEM
:
    public environmentMixingModel
{
    // Private data

        //- Turbulence model
        const compressible::turbulenceModel& flTurb_;

        //- Turbulent kinetic energy field
        const volScalarField& k_;

        //- Turbulent dissipation rate field
        const volScalarField& epsilon_;


public:

        //- Runtime type information
        TypeName("IEM");


    // Constructors

        //- Construct from components
        IEM
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~IEM();


    // Member Functions

        //- Mixing kernel
        virtual tmp<fvScalarMatrix> K
        (
            const volScalarField& meanMoment,
            const volScalarField& meanMomentVariance,
            const volScalarField& meanMixtureFraction
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace environmentMixingModel
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
