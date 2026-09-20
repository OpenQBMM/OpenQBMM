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
    Foam::PDFTransportModels::mixingModels::conditionalMixing

Description
    Mixing and reaction of a mixture fraction and of the compositions
    conditioned on it.

    The distribution carries the moments of the mixture fraction and, for
    every composition, its mean at each node of the mixture fraction: the
    quadrature of conditional means. Both are on [0, 1], which the model
    gives the quadrature as the supports of its dimensions. The number of
    compositions is nCompositions when the dictionary gives it, and
    otherwise the number of species of the thermo of the case, which is
    what the chemistry reaction model needs it to be.

    With the temperature closure conditionalEnthalpy of the reaction, the
    distribution carries one coordinate more, after the compositions: the
    absolute enthalpy of the fluid, conditioned on the mixture fraction as
    the compositions are, on the whole real line. It is transported and
    mixed with them, and no reaction changes it, since the enthalpy of an
    environment is what a reaction at constant pressure conserves; its
    moments are therefore where the heat lost to a wall or to radiation
    would enter, which an enthalpy linear in the mixture fraction cannot
    carry. The quadrature must then be given one scalar more than the
    compositions.

    The moments are carried by the velocity of the fluid, since a node of
    the mixture fraction is a state of the fluid and not of anything moving
    through it, and mixed by the kernel of the case. The mixing of every
    moment is closed in the moments themselves, so the model needs no
    quadrature to compute it; IEM is written for the moments of any number
    of coordinates, which is what a case selects here.

    Where the flux that carries the moments is a mass flux, as it is in a
    solver of variable density, the distribution is mass-weighted: the
    weight of a node is the fraction of the mass of the fluid in its
    environment, the moments obey d (rho m)/d t + div(rho U m) = div(rho
    Gamma grad m) + rho S, and the means they give are mass-weighted
    means, which is what mass fractions and an enthalpy per unit mass
    average to. The mean density of the fluid is then the one whose
    specific volume is the mean of those of the environments, 1/rho =
    sum_i w_i/rho_i. The density is named by rho in the dictionary of the
    model, rho unless another name is given, and has to obey continuity
    with that flux when the moments are solved, which a solver obtains by
    solving its continuity equation first.

    The reaction is a model of its own, selected in the reaction
    sub-dictionary and none without it, and it acts on the nodes rather
    than on the moments: once the moments have been transported and
    mixed, and the quadrature updated, every node is an environment whose
    compositions the reaction model advances over the time step, after
    which the moments are rebuilt from the nodes. The rebuild is exact,
    since a reaction moves neither the weights nor the mixture fraction,
    so the step is the moment equations of the reaction integrated
    node by node, split from the mixing as the mixing is split from
    transport.

    \verbatim
    conditionalMixingCoeffs
    {
        nCompositions   5;              // or the species of the thermo
        mixingKernel    { mixingKernel IEM; Cphi 2; }
        diffusionModel  { mixingDiffusionModel none; }
        reaction        { reactionModel chemistry; ... }
        odeCoeffs       { ... }         // of the realizable ODE solver
    }
    \endverbatim

SourceFiles
    conditionalMixing.C

\*---------------------------------------------------------------------------*/

#ifndef conditionalMixing_H
#define conditionalMixing_H

#include "velocityPDFTransportModel.H"
#include "mixingModel.H"
#include "mixingKernel.H"
#include "mixingDiffusionModel.H"
#include "reactionModel.H"
#include "realizableOdeSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace mixingModels
{

/*---------------------------------------------------------------------------*\
                      Class conditionalMixing Declaration
\*---------------------------------------------------------------------------*/

class conditionalMixing
:
    public velocityPDFTransportModel,
    public mixingModel,
    private realizableOdeSolver<volVelocityMoment, volVelocityNode>
{
    typedef realizableOdeSolver<volVelocityMoment, volVelocityNode> odeType;

    // Private data

        //- Name of the distribution
        const word name_;

        //- The moment fields of the quadrature, as the kernels read them
        UPtrList<const volScalarField> momentFields_;

        //- Mixing kernel
        autoPtr<Foam::mixingSubModels::mixingKernel> mixingKernel_;

        //- Diffusion model
        autoPtr<Foam::mixingSubModels::mixingDiffusionModel> diffusionModel_;

        //- Reaction model
        autoPtr<Foam::mixingSubModels::reactionModel> reactionModel_;


    // Private member functions

        //- Update implicit source terms
        virtual void updateImplicitMomentSource();

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volVelocityMoment& moment
        );

        //- Calculate explicit source terms with realizable ODE solver
        virtual void explicitMomentSource();

        //- Are moments updated due to sources
        virtual bool solveMomentSources() const;

        //- Is the realizable ode solver used to solve sources
        virtual bool solveMomentOde() const;

        //- Update explicit source terms
        virtual void updateCellMomentSource(const label celli);

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const labelList& momentOrder,
            const label celli,
            const velocityQuadratureApproximation& quadrature,
            const label environment
        );


public:

    //- Runtime type information
    TypeName("conditionalMixing");


    // Constructors

        //- Construct from components
        conditionalMixing
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~conditionalMixing();


    // Member Functions

        //- The number of compositions conditioned on the mixture fraction:
        //  nCompositions of the dictionary, or the species of the thermo
        //  of the case when the dictionary does not give it
        static label nCompositions
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Whether the distribution carries the enthalpy conditioned on
        //  the mixture fraction: the temperature closure of the reaction
        //  is conditionalEnthalpy
        static bool carriesEnthalpy(const dictionary& dict);

        //- The supports of the dimensions of the distribution: the mixture
        //  fraction and every composition on [0, 1], and the enthalpy,
        //  when it is carried, on the real line
        static List<supportType> supports
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- The reaction model
        const Foam::mixingSubModels::reactionModel& reaction() const
        {
            return reactionModel_();
        }

        //- Return the maximum Courant number ensuring moment realizability
        scalar realizableCo();

        //- Transport and mix the moments, then react the nodes
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingModels
} // End namespace PDFTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
