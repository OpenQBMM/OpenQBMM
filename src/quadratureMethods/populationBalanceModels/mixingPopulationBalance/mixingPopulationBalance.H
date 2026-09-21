/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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
    Foam::populationBalanceModels::mixingPopulationBalance

Description
    Population balance of a fluid that is mixing, in the two environments
    of a mixture fraction: the two nodes of the quadrature of a
    turbulentMixing model, of weights p_1, p_2 and mixture fractions xi_1,
    xi_2, each of which holds a population of its own, m_1 and m_2.

    Two sets of moments of the population are transported, the mean ones
    and the ones weighted by the mixture fraction (the quadratures named
    MeanMoments and MeanMomentsVariance),

    \f[
        \langle m \rangle = p_1 m_1 + p_2 m_2, \quad
        \langle \xi m \rangle = p_1 \xi_1 m_1 + p_2 \xi_2 m_2
    \f]

    the second with the source of the environmentMixingModel. At every step
    the populations of the environments are found from the two, changed by
    the sources of the population balance in each environment, with the
    realizable ODE solver, and put together again. This is done with the
    environments of the time the step starts from, which is the time of the
    moments: the mixing model is advanced after it.

    The environments are told apart where the variance of the mixture
    fraction, p_1 p_2 (xi_1 - xi_2)^2, is above minMixtureFractionVariance
    (1e-4) and each holds more than minEnvironmentWeight (1e-6) of the
    fluid. Elsewhere the fluid of the cell is mixed to what the thresholds
    resolve, both environments hold the mean moments, and <xi m> becomes
    <xi><m>: the covariance of such a cell is given up, so that the
    thresholds are to be small against the variances of interest. A cell
    where <xi m>/<m> lies outside the two mixture fractions, which would
    give an environment a negative population, is taken as mixed as well.

SourceFiles
    mixingPopulationBalance.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPopulationBalance_H
#define mixingPopulationBalance_H

#include "univariatePDFTransportModel.H"
#include "populationBalanceModel.H"
#include "turbulentMixing.H"
#include "aggregationKernel.H"
#include "breakupKernel.H"
#include "daughterDistribution.H"
#include "growthModel.H"
#include "diffusionModel.H"
#include "nucleationModel.H"
#include "environmentMixingModel.H"
#include "realizableOdeSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{

/*---------------------------------------------------------------------------*\
                    Class mixingPopulationBalance Declaration
\*---------------------------------------------------------------------------*/

class mixingPopulationBalance
:
    public populationBalanceModel,
    private realizableOdeSolver<volScalarMoment, volScalarNode>
{
    typedef realizableOdeSolver<volScalarMoment, volScalarNode> odeType;

    // Private data

        //- Turbulent mixing model
        autoPtr
        <
            Foam::PDFTransportModels::mixingModels::turbulentMixing
        > mixingModel_;

        //- Minimum mixture fraction variance to solve multienvironment problem
        scalar minMixtureFractionVariance_;

        //- Minimum environment weight
        scalar minEnvironmentWeight_;

        //- Reference to field of weight of environment one
        const volScalarField& p1_;

        //- Reference to field of abscissa of environment one
        const volScalarField& xi1_;

        //- Reference to field of weight of environment two
        const volScalarField& p2_;

        //- Reference to field of abscissa of environment two
        const volScalarField& xi2_;

        //- Reference to mean mixture fraction
        const volScalarField& meanXi_;

        //- Quadrature for the mean moments
        scalarQuadratureApproximation meanMomentsQuadrature_;

        //- Quadrature for the mean moment variance
        scalarQuadratureApproximation meanMomentsVarianceQuadrature_;

        //- Moment advection for mean moments
        autoPtr<univariateMomentAdvection> meanMomentsAdvection_;

        //- Moment advection for variance of mean moments
        autoPtr<univariateMomentAdvection> meanMomentsVarianceAdvection_;

        //- Mean moments
        volScalarMomentFieldSet& meanMoments_;

        //- Variance of mean moments
        volScalarMomentFieldSet& meanMomentsVariance_;

        //- Quadrature approximation of the PBE in the first environment
        scalarQuadratureApproximation envOneQuadrature_;

        //- Quadrature approximation of the PBE in the second environment
        scalarQuadratureApproximation envTwoQuadrature_;

        //- Reference to moments in environment one
        volScalarMomentFieldSet& mEnvOne_;

        //- Reference to moments in environment two
        volScalarMomentFieldSet& mEnvTwo_;

        //- Switch indicating the status of the aggregation model
        Switch aggregation_;

        //- Switch indicating the status of the breakup model
        Switch breakup_;

        //- Switch indicating the status of the growth model
        Switch growth_;

        //- Switch indicating the status of the nucleation model
        Switch nucleation_;

        //- Aggregation kernel model
        autoPtr
        <
            Foam::populationBalanceSubModels::aggregationKernel
        > aggregationKernel_;

        //- Breakup kernel model
        autoPtr
        <
            Foam::populationBalanceSubModels::breakupKernel
        > breakupKernel_;

        //- Growth model function
        autoPtr
        <
            Foam::populationBalanceSubModels::growthModel
        > growthModel_;

        //- Diffusion model
        autoPtr
        <
            Foam::populationBalanceSubModels::diffusionModel
        > diffusionModel_;

        //- Nucleation model
        autoPtr
        <
            Foam::populationBalanceSubModels::nucleationModel
        > nucleationModel_;

        //- Environment mixing model
        autoPtr
        <
            Foam::populationBalanceSubModels::environmentMixingModel
        > envMixingModel_;


    // Private member functions

        //- Calculate moments of the PBE in each environment
        void calcEnvironmentMoments();

        //- Calculated mixed moments
        void calcMixedMoments();

        //- Update explicit collison source terms
        virtual void updateCellMomentSource(const label celli);

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const labelList& momentOrder,
            const label celli,
            const scalarQuadratureApproximation& nodes,
            const label enviroment
        );


public:

    //- Runtime type information
    TypeName("mixingPbe");


    // Constructors

        //- Construct from components
        mixingPopulationBalance
        (
            const word& name,
            const dictionary& dict,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~mixingPopulationBalance();

    // Member Functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const;

        //- Return the courant number based on abscissae
        virtual scalar CoNum() const;

        //- Solve population balance equation
        void solve();

        //- Read
        virtual bool readIfModified();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace univariatePDFTransportModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
