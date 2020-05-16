/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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
    Foam::mixingSubModels::mixingDiffusionModels::turbulentDiffusion

Description
    Turbulent diffusion model.

    The moment diffusion is modeled as

    \f[
        \mathcal{D}_k = \nabla \cdot \left[ \left( \Gamma_l
            + \frac{\mu_t}{Sc_t} \right) \nabla m_k \right]
    \f]

    where

    \vartable
        \mathcal{D}_k       |   Diffusion term in the moment transport equation
        m_k                 |   Moment of order k [-]
        Sc_t                |   Turbulent Schmidt number [-]
        \Gamma_l            |   Molecular diffusivity [m2/s]
        \mu_t               |   Turbulent dynamic viscosity [m2/s]
    \endvartable

SourceFiles
    turbulentDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDiffusion_H
#define turbulentDiffusion_H

#include "mixingDiffusionModel.H"
#include "fvm.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace mixingDiffusionModels
{

/*---------------------------------------------------------------------------*\
                    Class turbulentDiffusion Declaration
\*---------------------------------------------------------------------------*/

class turbulentDiffusion
:
    public mixingDiffusionModel
{
    //- Private member functions
        tmp<volScalarField> turbViscosity(const volScalarField& moment) const;

    //- Private data

        //- Laminar diffusivity
        dimensionedScalar gammaLam_;

        //- Turbulent Schmidt number
        scalar Sc_;

public:

        //- Runtime type information
        TypeName("turbulentDiffusion");


    // Constructors

        //- Construct from components
        turbulentDiffusion(const dictionary& dict);


    //- Destructor
    virtual ~turbulentDiffusion();


    // Member Functions

        //- Moments of the daughter distribution
        virtual tmp<fvScalarMatrix> momentDiff
        (
            const volScalarField& moment
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingDiffusionModels
} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
