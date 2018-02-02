/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
     \\/     M anipulation  |
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
    Foam::populationBalanceSubModels::aggregationKernels::ChengVigilFoxFNP

Description
    Kernel for flash nano-precipitation.

    References
    \verbatim
        "A competitive aggregation model for flash nanoprecipitation"
        J C Cheng, R D Vigil, R O Fox
        Journal of Colloid and Interface Science
        Volume 351, Pages 330-342, 2010
    \endverbatim

SourceFiles
    ChengVigilFoxFNP.C

\*---------------------------------------------------------------------------*/

#ifndef ChengVigilFoxFNP_H
#define ChengVigilFoxFNP_H

#include "aggregationKernel.H"
#include "turbulentFluidThermoModel.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{

/*---------------------------------------------------------------------------*\
                    Class ChengVigilFoxFNP Declaration
\*---------------------------------------------------------------------------*/

class ChengVigilFoxFNP
:
    public aggregationKernel
{
    // Private data

        //- Thermo model
        const fluidThermo& flThermo_;

        //- Turbulence model
        const compressible::turbulenceModel& flTurb_;

        //- Temperature field
        const volScalarField& T_;

        //- Density field
        const volScalarField& rho_;

        //- Dynamic viscosity field
        const volScalarField& mu_;

        //- Turbulent dissipation rate field
        const volScalarField& epsilon_;

        //- Quadrature approximation of the mixing model
        const univariateQuadratureApproximation& mixingQuadrature_;

        //- Mixture fraction field
        const volScalarField& mixtureFraction_;

        //- Molar mass
        dimensionedScalar molarMass_;

        //- Density of species 1
        dimensionedScalar molarVol1_;

        //- Density of species 2
        dimensionedScalar molarVol2_;

        //- Initial solute concentration
        dimensionedScalar soluteConc_;

        //- Aggregation efficiency
        dimensionedScalar aggregationEfficiency_;

public:

        //- Runtime type information
        TypeName("ChengVigilFoxFNP");


    // Constructors

        //- Construct from components
        ChengVigilFoxFNP
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~ChengVigilFoxFNP();


    // Member Functions

        //- Aggregation kernel
        virtual scalar Ka
        (
            const scalar& abscissa1,
            const scalar& abscissa2,
            const label celli
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
