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
    Foam::populationBalanceSubModels::environmentMixingModel::IEM

Description
    Interaction by exchange with the mean mixing model, for the moments of
    a population weighted by the mixture fraction, <xi m>, in a population
    balance of the environments of a mixture fraction.

    The mixture fraction of an environment and the population it holds are
    drawn to their means at the same rate, Cphi/2 epsilon/k, which gives
    the source of <xi m>:

    \f[
        \frac{d \langle \xi m \rangle}{dt} =
        C_\phi \frac{\varepsilon}{k}
        \left(
            \langle \xi \rangle \langle m \rangle
          - \langle \xi m \rangle
        \right)
    \f]

    that is, the covariance of the population and the mixture fraction
    decays at the rate of the variance of the mixture fraction, and a
    population that is the mixture fraction itself keeps <xi m> at the
    second moment of it.

    References
    \verbatim
        "Computational Models for Turbulent Reacting Flows"
        R O Fox
        Cambridge University Press
    \endverbatim

SourceFiles
    IEM.C

\*---------------------------------------------------------------------------*/

#ifndef IEM_H
#define IEM_H

#include "environmentMixingModel.H"
#include "turbulentFluidThermoModel.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace environmentMixingModels
{

/*---------------------------------------------------------------------------*\
                           Class IEM Declaration
\*---------------------------------------------------------------------------*/

class IEM
:
    public environmentMixingModel
{
    // Private data

        //- Turbulence model
        const compressible::turbulenceModel& flTurb_;

        //- Turbulent kinetic energy and its dissipation rate, as the
        //  turbulence model returned them at the time step they were last
        //  asked for. A model that does not solve for one of them builds it
        //  anew each time, so they are held, not referred to.
        mutable tmp<volScalarField> k_;
        mutable tmp<volScalarField> epsilon_;

        //- Time step k_ and epsilon_ were taken at
        mutable label turbulenceTimeIndex_;


    // Private Member Functions

        //- Take k and epsilon from the turbulence model, once per time step
        void updateTurbulence() const;

        //- Turbulent kinetic energy
        const volScalarField& k() const;

        //- Turbulent dissipation rate
        const volScalarField& epsilon() const;


public:

        //- Runtime type information
        TypeName("IEM");


    // Constructors

        //- Construct from components
        IEM
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~IEM();


    // Member Functions

        //- Mixing kernel
        virtual tmp<fvScalarMatrix> K
        (
            const volScalarField& meanMoment,
            const volScalarField& meanMomentVariance,
            const volScalarField& meanMixtureFraction
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace environmentMixingModel
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
