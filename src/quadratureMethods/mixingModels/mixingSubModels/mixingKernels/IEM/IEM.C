/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 Alberto Passalacqua
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

#include "IEM.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace mixingKernels
{
    defineTypeNameAndDebug(IEM, 0);

    addToRunTimeSelectionTable
    (
        mixingKernel,
        IEM,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::IEM::IEM
(
    const dictionary& dict
)
:
    mixingKernel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::IEM
::~IEM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::mixingSubModels::mixingKernels::IEM::K
(
    const volMoment& moment,
    const volMomentFieldSet& moments
) const
{
    typedef compressible::turbulenceModel cmpTurbModel;

    if
    (
        !moment.mesh().foundObject<cmpTurbModel>
        (
            cmpTurbModel::propertiesName
        )
    )
    {
        FatalErrorInFunction
            << "No valid compressible turbulence model found."
            << abort(FatalError);
    }

    const compressible::turbulenceModel& flTurb =
        moment.mesh().lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    label momentOrder = moment.order();

    tmp<fvScalarMatrix> mixingK
    (
        new fvScalarMatrix
        (
            moment,
            moment.dimensions()*dimVol/dimTime
        )
    );

    if (momentOrder == 0)
    {
        return mixingK;
    }
    else
    {
        mixingK.ref() += momentOrder*Cphi_*flTurb.epsilon()/flTurb.k()
            *(moments[momentOrder - 1]*moments[1])
            - fvm::SuSp(momentOrder*Cphi_*flTurb.epsilon()/flTurb.k(), moment);
    }

    return mixingK;
}

// ************************************************************************* //
