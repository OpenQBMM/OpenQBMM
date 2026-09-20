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
    Foam::mixingSubModels::reactionModels::noReaction

Description
    No reaction: the compositions conditioned on the mixture fraction are
    only mixed. Selected as none, and what a mixing model uses when its
    dictionary has no reaction sub-dictionary.

SourceFiles
    noReaction.C

\*---------------------------------------------------------------------------*/

#ifndef noReaction_H
#define noReaction_H

#include "reactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace reactionModels
{

/*---------------------------------------------------------------------------*\
                         Class noReaction Declaration
\*---------------------------------------------------------------------------*/

class noReaction
:
    public reactionModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noReaction
        (
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            const label nCompositions
        );


    //- Destructor
    virtual ~noReaction();


    // Member Functions

        //- The compositions never change
        virtual bool reacts() const
        {
            return false;
        }

        //- Nothing to do
        virtual void react(velocityQuadratureApproximation& quadrature);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactionModels
} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
