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
    Foam::populationBalanceSubModels::growthModels::coolCrysSizeEffectGrowth

Description
    Crystal growth model with a size-effect inhibition factor.

SourceFiles
    coolCrysSizeEffectGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef coolCrysSizeEffectGrowth_H
#define coolCrysSizeEffectGrowth_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
            Class coolCrysSizeEffectGrowth Declaration
\*---------------------------------------------------------------------------*/

class coolCrysSizeEffectGrowth
:
    public growthModel
{
private:

    // Private data

        //- Growth rate coefficient
        dimensionedScalar Cg_;

        //- Power coefficient for growth law
        scalar powCoeff_;

        //- Critical supersaturation for growth onset
        scalar sigmaGrowthCrit_;

        //- Upper bound for supersaturation used in the growth law
        scalar sigmaMax_;

        //- Upper bound for the returned growth rate
        dimensionedScalar maxGrowth_;

        //- Characteristic size for the inhibition factor
        dimensionedScalar Lc_;

        //- Size-effect exponent
        scalar qSize_;

        //- Supersaturation field
        const volScalarField& sigma_;

        //- Growth rate field for output
        volScalarField GField_;


    // Private Member Functions

        //- Growth rate for a characteristic crystal length in a cell
        scalar growthRateForLength
        (
            const scalar& length,
            const label celli
        ) const;

        //- Update growthRate with a number-density weighted mean
        void updateMeanGrowthRate
        (
            const label celli,
            const scalarQuadratureApproximation& quadrature
        ) const;

        //- Update growthRate with a number-density weighted mean
        void updateMeanGrowthRate
        (
            const label celli,
            const velocityQuadratureApproximation& quadrature
        ) const;


public:

        //- Runtime type information
        TypeName("coolCrysSizeEffectGrowth");


    // Constructors

        //- Construct from components
        coolCrysSizeEffectGrowth
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~coolCrysSizeEffectGrowth();


    // Member Functions

        //- Growth model
        virtual scalar Kg
        (
            const scalar& abscissa,
            const bool lengthBased,
            const label celli,
            const label environment = 0
        ) const;

        //- Return phase space convection and update mean growthRate output
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& quadrature
        );

        //- Return phase space convection and update mean growthRate output
        virtual scalar phaseSpaceConvection
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
