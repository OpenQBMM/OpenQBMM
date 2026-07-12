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
    Foam::populationBalanceSubModels::growthModels::coolCrysAreaInhibitionGrowth

Description
    Crystal growth model with inhibition driven by the second-moment field.

SourceFiles
    coolCrysAreaInhibitionGrowth.C

\*---------------------------------------------------------------------------*/

#ifndef coolCrysAreaInhibitionGrowth_H
#define coolCrysAreaInhibitionGrowth_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace growthModels
{

/*---------------------------------------------------------------------------*\
            Class coolCrysAreaInhibitionGrowth Declaration
\*---------------------------------------------------------------------------*/

class coolCrysAreaInhibitionGrowth
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

        //- Name of the second-moment field
        word m2Name_;

        //- Area inhibition coefficient
        dimensionedScalar alphaArea_;

        //- Exponent in the area-inhibition factor
        scalar areaPower_;

        //- Supersaturation field
        const volScalarField& sigma_;

        //- Cell being evaluated by phaseSpaceConvection
        label activeCelli_;

        //- Second-moment field used for inhibition
        const volScalarField* m2Ptr_;

        //- Growth rate field for output
        volScalarField GField_;


public:

        //- Runtime type information
        TypeName("coolCrysAreaInhibitionGrowth");


    // Constructors

        //- Construct from components
        coolCrysAreaInhibitionGrowth
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~coolCrysAreaInhibitionGrowth();


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
