/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015 by Matteo Icardi
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::populationBalanceSubModels::diffusionModels::molecularDiffusion

Description
    Molecular diffusion model.

    The moment diffusion is modeled with a constant diffusion coefficient
    specified by the user.

    \f[
        \mathcal{D}_k = \nabla \cdot \left( \Gamma_l
            \nabla m_k \right)
    \f]

    where

    \vartable
        \mathcal{D}_k       |   Diffusion term in the moment transport equation
        m_k                 |   Moment of order k [-]
        \Gamma_l            |   Molecular diffusivity [m2/s]
    \endvartable

SourceFiles
    molecularDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef molecularDiffusion_H
#define molecularDiffusion_H

#include "diffusionModel.H"
#include "fvm.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace diffusionModels
{

/*---------------------------------------------------------------------------*\
                    Class molecularDiffusion Declaration
\*---------------------------------------------------------------------------*/

class molecularDiffusion
:
    public diffusionModel
{
    //- Private member functions
        tmp<volScalarField> turbViscosity(const volScalarField& moment) const;

    //- Private data

        //- Laminar diffusivity
        dimensionedScalar gammaLam_;

public:

        //- Runtime type information
        TypeName("molecularDiffusion");


    // Constructors

        //- Construct from components
        molecularDiffusion(const dictionary& dict);


    //- Destructor
    virtual ~molecularDiffusion();


    // Member Functions

        //- Moments of the daughter distribution
        virtual tmp<fvScalarMatrix> momentDiff
        (
            const volScalarField& moment
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diffusionModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
