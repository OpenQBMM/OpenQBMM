/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::nucleationModels::noNucleation

Description
    No nucleation model.

SourceFiles
    noNucleation.C

\*---------------------------------------------------------------------------*/

#ifndef noNucleation_H
#define noNucleation_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace nucleationModels
{

/*---------------------------------------------------------------------------*\
                        Class noNucleation Declaration
\*---------------------------------------------------------------------------*/

class noNucleation
:
    public nucleationModel
{
public:

        //- Runtime type information
        TypeName("none");


    // Constructors

        //- Construct from components
        noNucleation
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~noNucleation();


    // Member Functions

        //- Growth model
        virtual tmp<volScalarField> nucleationSource
        (
            const volUnivariateMoment& moment
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationModels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
