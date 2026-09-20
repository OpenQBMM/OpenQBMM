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
    const dictionary& dict,
    const fvMesh& mesh,
    const labelListList& momentOrders,
    const UPtrList<const volScalarField>& moments
)
:
    mixingKernel(dict, mesh, momentOrders, moments)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::IEM::~IEM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::mixingSubModels::mixingKernels::IEM::mixingSource
(
    const labelList& momentOrder,
    const label celli,
    const label environment
) const
{
    // Every value relaxes to its mean at the rate Cphi/2 epsilon/k, so the
    // variance decays at Cphi epsilon/k. The half was missing, which with
    // the usual Cphi = 2 mixed twice as fast as the model does.
    //
    // The moment of orders k of any number of coordinates then has the
    // source -w sum_i k_i (m(k) - <phi_i> m(k - e_i)), in which every term
    // is a moment the distribution carries: the mixing of a mixture
    // fraction and of the compositions conditioned on it is closed in its
    // moments, and needs no quadrature.
    const scalar m0 = moment(labelList(momentOrder.size(), 0), celli);

    if (m0 <= 0)
    {
        return 0.0;
    }

    const scalar omega = 0.5*Cphi_.value()*epsilon()[celli]/k()[celli];
    const scalar m = moment(momentOrder, celli);

    scalar source = 0;

    forAll(momentOrder, dimi)
    {
        if (momentOrder[dimi] == 0)
        {
            continue;
        }

        labelList firstOrder(momentOrder.size(), 0);
        firstOrder[dimi] = 1;

        const scalar mean = moment(firstOrder, celli)/m0;

        source +=
            momentOrder[dimi]
           *(mean*momentBelow(momentOrder, dimi, celli) - m);
    }

    return omega*source;
}

// ************************************************************************* //
