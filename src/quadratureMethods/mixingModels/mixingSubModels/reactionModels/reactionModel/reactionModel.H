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
    Foam::mixingSubModels::reactionModel

Description
    Base class of the reaction of the compositions conditioned on the
    mixture fraction, one environment at a time.

    The quadrature of conditional means carries n nodes of the mixture
    fraction, each with the mean of every composition conditioned on it:
    an environment is a node, a weight w_i, a mixture fraction xi_i and a
    composition Y_i. A reaction changes the compositions of an environment
    and nothing else: no reaction moves the mixture fraction, which is
    conserved, or the weight of a node. The source of the moment of orders
    (k, e_j), the moment of xi^k and of the composition j, is therefore

        d<xi^k Y_j>/dt = sum_i w_i xi_i^k omega_j(Y_i, T_i, p),

    the rate of the composition j in every environment, weighted by the
    node, which the quadrature gives without any closure.

    A reaction model acts on the nodes rather than on the moments: it
    advances the composition of every environment over the time step, with
    the weights and the mixture fraction as they are, and the moments are
    rebuilt from the nodes. The rebuild is exact, since the moments are
    sums over the nodes of products of quantities the reaction leaves
    untouched with the compositions it advanced; integrating the nodes and
    integrating the moment equations above are the same thing. The
    reaction is a separate step after the mixing of the moments, which is
    the same first-order splitting the laminar combustion model of
    OpenFOAM applies between transport and chemistry.

    What every environment reacts with is the model's: the concrete model
    says what compositions are, how many the quadrature must carry, how the
    temperature of an environment is found, and what it does with the
    mean composition of a cell once the environments have reacted.

    Selected in the reaction sub-dictionary of the mixing model by
    reactionModel, none when the sub-dictionary is missing:

    \verbatim
    reaction
    {
        reactionModel   chemistry;   // none | chemistry
        ...                          // entries of the model
    }
    \endverbatim

SourceFiles
    reactionModel.C
    newReactionModel.C

\*---------------------------------------------------------------------------*/

#ifndef reactionModel_H
#define reactionModel_H

#include "dictionary.H"
#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{

/*---------------------------------------------------------------------------*\
                       Class reactionModel Declaration
\*---------------------------------------------------------------------------*/

class reactionModel
{
    // Private member functions

        //- No copy construct
        reactionModel(const reactionModel&) = delete;

        //- No copy assignment
        void operator=(const reactionModel&) = delete;


protected:

    // Protected data

        //- Dictionary of the model
        const dictionary dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Number of compositions conditioned on the mixture fraction
        const label nCompositions_;


public:

    //- Runtime type information
    TypeName("reactionModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        reactionModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const label nCompositions
        ),
        (dict, mesh, nCompositions)
    );


    // Constructors

        //- Construct from the dictionary of the reaction, the mesh and the
        //  number of compositions the quadrature carries
        reactionModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const label nCompositions
        );


    // Selectors

        //- Select the model named by reactionModel in the dictionary, none
        //  when the dictionary has no such entry
        static autoPtr<reactionModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const label nCompositions
        );


    //- Destructor
    virtual ~reactionModel();


    // Member Functions

        //- Number of compositions conditioned on the mixture fraction
        inline label nCompositions() const
        {
            return nCompositions_;
        }

        //- Whether the model changes the compositions
        virtual bool reacts() const = 0;

        //- Advance the compositions of every environment of every cell
        //  over the time step of the case, and rebuild the moments from
        //  the nodes
        virtual void react(velocityQuadratureApproximation& quadrature) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
