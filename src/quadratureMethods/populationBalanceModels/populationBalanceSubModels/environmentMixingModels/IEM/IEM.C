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

\*---------------------------------------------------------------------------*/

#include "IEM.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace environmentMixingModels
{
    defineTypeNameAndDebug(IEM, 0);

    addToRunTimeSelectionTable
    (
        environmentMixingModel,
        IEM,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::environmentMixingModels::IEM::IEM
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    environmentMixingModel(dict, mesh),
    flTurb_
    (
        mesh_.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    ),
    k_(),
    epsilon_(),
    turbulenceTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::environmentMixingModels::IEM::~IEM()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
Foam::populationBalanceSubModels::environmentMixingModels::IEM
::updateTurbulence() const
{
    // This used to bind a reference to what the turbulence model returned
    // on construction. A model that solves for the field returns its own,
    // which the reference followed; one that does not, such as the k-omega
    // family for epsilon, returns a field built for the call, which the
    // reference outlived.
    const label timeIndex = mesh_.time().timeIndex();

    if (timeIndex != turbulenceTimeIndex_)
    {
        k_ = flTurb_.k();
        epsilon_ = flTurb_.epsilon();
        turbulenceTimeIndex_ = timeIndex;
    }
}


const Foam::volScalarField&
Foam::populationBalanceSubModels::environmentMixingModels::IEM::k() const
{
    updateTurbulence();
    return k_();
}


const Foam::volScalarField&
Foam::populationBalanceSubModels::environmentMixingModels::IEM::epsilon() const
{
    updateTurbulence();
    return epsilon_();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::populationBalanceSubModels::environmentMixingModels::IEM::K
(
    const volScalarField& meanMoment,
    const volScalarField& meanMomentVariance,
    const volScalarField& meanMixtureFraction
) const
{
    // Unverified, left as it was: the interaction by exchange with the mean
    // acting on the mixture fraction relaxes <xi m> to <xi><m> at
    // Cphi/2 epsilon/k, which is a quarter of the rate written here. The
    // mixing kernels had lost a factor of one half, now restored; whether
    // this model carries a reason of its own for the rest is to be checked
    // against its derivation before the two-environment model is used.
    return
        2.0*Cphi_*epsilon()*meanMoment*meanMixtureFraction/k()
      - fvm::SuSp(2.0*Cphi_*epsilon()/k(), meanMomentVariance);
}

// ************************************************************************* //
