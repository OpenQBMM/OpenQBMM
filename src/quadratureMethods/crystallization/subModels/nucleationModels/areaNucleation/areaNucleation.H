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
    Copyright (C) 2019-2023 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::nucleationModels::areaNucleation

Description
    Nucleation model with surface-area inhibition.

SourceFiles
    areaNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef areaNucleation_H
#define areaNucleation_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                        Class areaNucleation Declaration
\*---------------------------------------------------------------------------*/

class areaNucleation
:
    public nucleationModel
{
private:

    // Private data

        //- Nucleation rate coefficient k_b
        dimensionedScalar kb_;

        //- Power law exponent for supersaturation
        scalar b_;

        //- Surface-area inhibition coefficient
        scalar alpha_;

        //- Group name used for moment lookup
        word momentGroup_;

        //- Supersaturation field
        const volScalarField& sigma_;

        //- Second moment field
        const volScalarField& m2_;

        //- Nucleation rate field for output
        volScalarField JField_;

        //- Nucleation size
        dimensionedScalar d_nucleation_;


public:

    //- Runtime type information
    TypeName("areaNucleation");


    // Constructors

        //- Construct from components
        areaNucleation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~areaNucleation();


    // Member Functions

        //- Nucleation source term
        virtual scalar nucleationSource
        (
            const label& momentOrder,
            const label celli,
            const label environment = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
