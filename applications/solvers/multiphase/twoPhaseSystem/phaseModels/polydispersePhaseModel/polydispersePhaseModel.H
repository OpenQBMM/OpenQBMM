/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::polydispersePhaseModel

Description
    Class for a polydisperse phase model. Each bubble size hase a unique mean
    velocity. Size and velocity moments are stored in quadrature.

SourceFiles
    polydispersePhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef polydispersePhaseModel_H
#define polydispersePhaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "phaseModel.H"
#include "monoKineticQuadratureApproximation.H"

#include "coalescenceKernel.H"
#include "breakupKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class polydispersePhaseModel
:
    public phaseModel
{
    // Private data

        //- Population balance properties
        IOdictionary pbeDict_;

        //- Switch for ode source
        Switch solveOde_;

        //- Switch for coalescence
        Switch coalescence_;

        //- Switch for breakup
        Switch breakup_;

        //- QuadratureApproximation
        monoKineticQuadratureApproximation quadrature_;

        //- Number of nodes
        const label nNodes_;

        //- Number of moments
        const label nMoments_;

        //- List of alphas
        PtrList<volScalarField> alphas_;

        //- List of velocities (Refrence to the velocities stored in the
        //  quadrature
        PtrList<volVectorField>& Us_;

        //- List of deviation velocities
        PtrList<volVectorField> Vs_;

        //- Scale factor used for diameter
        scalar diameterScale_;

        //- List of diameters
        PtrList<volScalarField> ds_;

        //- Maximum diameter
        dimensionedScalar maxD_;

        //- Minimum diameter
        dimensionedScalar minD_;

        //- Correction field used to keep m1 and alpha*rho consistent
        tmp<volScalarField> corr_;

        //- Coalescence kernel
        autoPtr
        <
            populationBalanceSubModels::aggregationKernels::coalescence
        > coalescenceKernel_;

        //- Breakup kernel
        autoPtr<populationBalanceSubModels::breakupKernel> breakupKernel_;


        // ODE Coefficients

            //- Minimum local ode delta t
            scalar minLocalDt_;

            //- Stored field of time steps
            scalarField localDt_;

            //- Absolute tolerance of the abscissae
            scalar ATol_;

            //- Relative tolerance
            scalar RTol_;

            //- maximum increase factor in time step
            scalar facMax_;

            //- Maximum decrease factor in time step
            scalar facMin_;

            //- Precentage of error based time step
            scalar fac_;

            //- Vector to remove buildup due to coalescence and breakup
            //  in unused directions
            vector validDirections_;


    // Private member functions

        //- Update moments
        void updateVelocity();

        //- Coalesence source
        scalar coalescenceSource(const label celli, const label momentOrder);

        //- Coalesence source
        vector coalescenceSourceU(const label celli, const label momentOrder);

        //- Breakup source
        scalar breakupSource(const label celli, const label momentOrder);

        //- Breakup source
        vector breakupSourceU(const label celli, const label momentOrder);

        //- Solve ode for breakup and coalescnce
        void solveSourceOde();


public:

    // Constructors
        polydispersePhaseModel
        (
            const twoPhaseSystem& fluid,
            const dictionary& phaseProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~polydispersePhaseModel();


    // Member Functions

        //- Set population balance models
        virtual void setModels();

        //- Return the number of nodes
        virtual label nNodes() const
        {
            return nNodes_;
        }

        //- Return alpha field for nodei
        virtual const volScalarField& alphas(const label nodei) const
        {
            return alphas_[nodei];
        }

        //- Return non-constant access to alpha field for nodei
        virtual volScalarField& alphas(const label nodei)
        {
            return alphas_[nodei];
        }

        //- Return the diameter for nodei
        virtual const volScalarField& ds(const label nodei) const
        {
            return ds_[nodei];
        }

        //- Return the velocity for nodei
        virtual const volVectorField& Us(const label nodei) const
        {
            return Us_[nodei];
        }

        //- Return non-const access to the velocity for nodei
        virtual volVectorField& Us(const label nodei)
        {
            return Us_[nodei];
        }

        //- Return deviation velocity for nodei
        virtual tmp<volVectorField> Vs(const label nodei) const
        {
            return Vs_[nodei];
        }

        //- Correct the phase properties
        virtual void correct();

        //- Relative transport of moments
        virtual void relativeTransport();

        //- Average transport of moments
        virtual void averageTransport(const PtrList<fvVectorMatrix>& AEqns);

        //- Read phase properties dictionary
        virtual bool read(const bool readOK);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
