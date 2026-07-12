/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Junjie Li
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
    Foam::populationBalanceSubModels::aggregationKernels::
    crystalCollisionFrequencies::turbulentFrequency

Description
    Turbulence-induced collision frequency model for crystal aggregation.

    The collision frequency is calculated as:
    \f[
        \beta_{turb}(L_i, L_j) = C_{turb} \left(\frac{\varepsilon}{\nu}\right)^{1/2}
            (L_i + L_j)^3
    \f]

    \vartable
        C_{turb}    | Turbulent collision coefficient
        \varepsilon | Turbulent dissipation rate [m^2/s^3]
        \nu         | Kinematic viscosity [m^2/s]
        L_i, L_j    | Crystal sizes [m]
    \endvartable

SourceFiles
    turbulentFrequency.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentFrequency_H
#define turbulentFrequency_H

#include "crystalCollisionFrequency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
namespace crystalCollisionFrequencies
{

/*---------------------------------------------------------------------------*\
                    Class turbulentFrequency Declaration
\*---------------------------------------------------------------------------*/

class turbulentFrequency
:
    public crystalCollisionFrequency
{
    // Private data

        //- Turbulent collision coefficient
        //  Default value based on Saffman-Turner model: 1.91 * alpha
        //  where alpha ~ 0.797 for isotropic turbulence
        scalar Cturb_;

        //- Reference to kinematic viscosity field
        const volScalarField& nu_;

        //- Reference to turbulent dissipation rate
        const volScalarField& epsilon_;


public:

    //- Runtime type information
    TypeName("turbulent");


    // Constructors

        //- Construct from components
        turbulentFrequency
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~turbulentFrequency();


    // Member Functions

        //- Calculate collision frequency
        virtual scalar beta
        (
            const scalar& d1,
            const scalar& d2,
            const vector& Ur,
            const label celli
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace crystalCollisionFrequencies
} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
