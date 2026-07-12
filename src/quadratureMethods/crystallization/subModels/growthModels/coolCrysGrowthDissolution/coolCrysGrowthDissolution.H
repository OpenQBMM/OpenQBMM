/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::growthModels::PruppacherKlett

Description
    Diffusion-driven growth model. The typical value for Cg is 0.78.
    
    Reference
    \verbatim
        "Microphysics of Clouds and Precipitation"
        H R Pruppacher, J D Klett
        Reidel, Boston, MA, 1980
    \endverbatim

SourceFiles
    coolCrysGrowthDissolution.C

\*---------------------------------------------------------------------------*/

#ifndef coolCrysGrowthDissolution_H
#define coolCrysGrowthDissolution_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
            Class coolCrysGrowthDissolution Declaration
\*---------------------------------------------------------------------------*/

class coolCrysGrowthDissolution
:
    public growthModel
{
private:

    // Private data


        //- Growth rate coefficient
        dimensionedScalar Cg_;

        //- Power coefficient for growth law
        scalar powCoeff_;

        dimensionedScalar CgDiss_;

        scalar powCoeffDiss_;
        
        //- Oversaturation field (driving force for growth)
        const volScalarField& sigma_;

        //- Cell being evaluated by phaseSpaceConvection
        label activeCelli_;
        
        //- Growth rate field for output
        volScalarField GField_;

    
public:

        //- Runtime type information
        TypeName("coolCrysGrowthDissolution");


    // Constructors

        //- Construct from components
        coolCrysGrowthDissolution
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~coolCrysGrowthDissolution();


    // Member Functions

        //- Growth model
        virtual scalar Kg
        (
            const scalar& abscissa,
            const bool lengthBased,
            const label environment = 0
        ) const override;

        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature
        ) override;

        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature
        ) override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
