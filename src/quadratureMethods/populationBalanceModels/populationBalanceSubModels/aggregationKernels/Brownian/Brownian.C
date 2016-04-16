/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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

#include "Brownian.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace aggregationKernels
{
    defineTypeNameAndDebug(Brownian, 0);

    addToRunTimeSelectionTable
    (
        aggregationKernel,
        Brownian,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::Brownian
::Brownian
(
    const dictionary& dict
)
:
    aggregationKernel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernels::Brownian
::~Brownian()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::aggregationKernels::Brownian::Ka
(
    const volScalarField& abscissa1,
    const volScalarField& abscissa2
) const
{
    if (!abscissa1.mesh().foundObject<fluidThermo>(basicThermo::dictName))
    {
        FatalErrorInFunction
            << "No valid thermophysical model found."
            << abort(FatalError);
    }

    const fluidThermo& flThermo =
        abscissa1.mesh().lookupObject<fluidThermo>(basicThermo::dictName);

    dimensionedScalar smallAbs("smallAbs", sqr(abscissa1.dimensions()), SMALL);

    return
        2.0*Foam::constant::physicoChemical::k*flThermo.T()
        *sqr(abscissa1 + abscissa2)/(3.0*flThermo.mu()
        *max(abscissa1*abscissa2, smallAbs));
}

// ************************************************************************* //
