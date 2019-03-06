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

#include "FokkerPlanck.H"
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
    defineTypeNameAndDebug(FokkerPlanck, 0);

    addToRunTimeSelectionTable
    (
        mixingKernel,
        FokkerPlanck,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::FokkerPlanck
::FokkerPlanck
(
    const dictionary& dict
)
:
    mixingKernel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::FokkerPlanck::~FokkerPlanck()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::mixingSubModels::mixingKernels::FokkerPlanck::K
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

    dimensionedScalar oneMoment("oneMoment", moments[1].dimensions(), 1.0);

    if (momentOrder == 0)
    {
        return mixingK;
    }
    else
    {
        mixingK.ref() += momentOrder*Cphi_*flTurb.epsilon()/flTurb.k()
            *moments[momentOrder - 1]
            *((Cmixing_ + 1.0)*moments[1] + Cmixing_*(momentOrder - 1)*oneMoment
            *((moments[2] - sqr(moments[1]))/(moments[1]*oneMoment
            - moments[2]))) - fvm::SuSp(momentOrder*Cphi_*flTurb.epsilon()
            /flTurb.k()*((Cmixing_ + 1.0) + Cmixing_*(momentOrder - 1)
            *((moments[2] - sqr(moments[1]))/(moments[1]*oneMoment
            - moments[2]))), moment);
    }

    return mixingK;
}

// ************************************************************************* //
