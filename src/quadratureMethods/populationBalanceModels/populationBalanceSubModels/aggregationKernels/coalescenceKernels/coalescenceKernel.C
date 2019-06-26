/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/

#include "coalescenceKernel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(coalescence, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        coalescence,
        dictionary
    );
}
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescence::
coalescence
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    aggregationKernel(dict, mesh),
    continuousPhase_(dict.lookupOrDefault("contiunousPhase", word::null)),
    frequency_
    (
        coalescenceFrequencyKernel::New(dict, mesh, continuousPhase_)
    ),
    efficiency_
    (
        coalescenceEfficiencyKernel::New(dict, mesh, continuousPhase_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::coalescence::
~coalescence()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::populationBalanceSubModels::aggregationKernels::coalescence::
preUpdate()
{
    const fluidThermo& thermo =
        mesh_.lookupObject<fluidThermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                continuousPhase_
            )
        );
    const compressible::turbulenceModel& turb =
        mesh_.lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                continuousPhase_
            )
        );

    frequency_->update(thermo, turb);
    efficiency_->update(thermo, turb);
}

Foam::scalar
Foam::populationBalanceSubModels::aggregationKernels::coalescence::Ka
(
    const scalar& d1,
    const scalar& d2,
    const vector& Ur,
    const label celli,
    const label environment
) const
{
    return
        Ca_.value()
       *frequency_->omega(d1, d2, Ur, celli)
       *efficiency_->Pc(d1, d2, Ur, celli);
}


// ************************************************************************* //
