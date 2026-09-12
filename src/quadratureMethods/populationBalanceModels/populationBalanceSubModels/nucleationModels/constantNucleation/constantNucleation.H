/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::nucleationModels::constantNucleation

Description
    Particles of one size, formed at a rate that does not change.

    The rate and the size are read from the dictionary, so the model suits
    any internal coordinate: the size is written in the units the abscissae
    of the case are in, a length, a volume or a mass.

SourceFiles
    constantNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef constantNucleation_H
#define constantNucleation_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                     Class constantNucleation Declaration
\*---------------------------------------------------------------------------*/

class constantNucleation
:
    public nucleationModel
{
    // Private data

        //- Number of particles formed per unit volume and time
        dimensionedScalar nucleationRate_;

        //- Size the particles are formed with
        dimensionedScalar nucleationSize_;


public:

        //- Runtime type information
        TypeName("constant");


    // Constructors

        //- Construct from components
        constantNucleation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~constantNucleation();


    // Member Functions

        //- Nucleation rate
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
