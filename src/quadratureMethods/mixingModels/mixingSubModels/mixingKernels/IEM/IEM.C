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
    Copyright (C) 2019-2020 Alberto Passalacqua
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
    const fvMesh& mesh
)
:
    mixingKernel(dict, mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::mixingKernels::IEM::~IEM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix>
Foam::mixingSubModels::mixingKernels::IEM::K
(
    const volScalarMoment& moment,
    const volScalarMomentFieldSet& moments
) const
{
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
        mixingK.ref() += momentOrder*Cphi_*epsilon()/k()
            *(moments(momentOrder - 1)*moments(1))
            - fvm::SuSp(momentOrder*Cphi_*epsilon()/k(), moment);
    }

    return mixingK;
}

// ************************************************************************* //
