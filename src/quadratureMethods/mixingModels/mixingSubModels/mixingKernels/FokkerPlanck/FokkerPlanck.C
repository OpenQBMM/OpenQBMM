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
    Copyright (C) 2019-2024 Alberto Passalacqua
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
    const dictionary& dict,
    const fvMesh& mesh,
    const labelListList& momentOrders,
    const UPtrList<const volScalarField>& moments
)
:
    mixingKernel(dict, mesh, momentOrders, moments)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::FokkerPlanck::~FokkerPlanck()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::mixingSubModels::mixingKernels::FokkerPlanck::mixingSource
(
    const labelList& momentOrder,
    const label celli,
    const label environment
) const
{
    forAll(momentOrder, dimi)
    {
        if (dimi > 0 && momentOrder[dimi] != 0)
        {
            FatalErrorInFunction
                << "The Fokker-Planck kernel is written for the moments of "
                << "one coordinate, and the moment of orders "
                << momentOrder << " is of several." << nl
                << "    Mix the coordinates conditioned on it with IEM."
                << nl << exit(FatalError);
        }
    }

    const label order = momentOrder[0];

    if (order == 0)
    {
        return 0.0;
    }
    
    // The drift relaxes each value to the mean at (1 + Cmixing) Cphi/2
    // epsilon/k and the diffusion puts back Cmixing Cphi epsilon/k of the
    // variance, so the variance decays at Cphi epsilon/k whatever Cmixing
    // is, as with IEM. The half was missing from the drift, which made the
    // decay (2 + Cmixing) Cphi epsilon/k.
    const scalar omega = Cphi_.value()*epsilon()[celli]/k()[celli];

    const scalar m1 = moment(labelList(1, 1), celli);
    const scalar m2 = moment(labelList(1, 2), celli);
    const scalar mk = moment(labelList(1, order), celli);
    const scalar mkBelow = moment(labelList(1, order - 1), celli);

    return
        omega*order
       *(
            0.5*(Cmixing_.value() + 1.0)*(mkBelow*m1 - mk)
          + (order - 1)*Cmixing_.value()*(m2 - sqr(m1))*(mkBelow - mk)
           /max(2.0*(m1 - m2), SMALL)
        );
}

// ************************************************************************* //
