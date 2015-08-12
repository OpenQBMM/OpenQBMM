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

#include "LuoSvendsen.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(LuoSvendsen, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        LuoSvendsen,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen
::LuoSvendsen
(
    const dictionary& dict
)
:
    breakupKernel(dict),
    Cb_(dict.lookup("Cb")),
    epsilonExp_(readScalar(dict.lookup("epsilonExp"))),
    nuExp_(readScalar(dict.lookup("nuExp"))),
    sizeExp_(readScalar(dict.lookup("sizeExp")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen::~LuoSvendsen()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::breakupKernels::LuoSvendsen::Kb
(
    const volScalarField& abscissa
) const
{   
    const compressible::turbulenceModel& fluidTurb =
        abscissa.mesh().lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    return Cb_*pow(fluidTurb.epsilon(), epsilonExp_)
        *pow(fluidTurb.nu(), nuExp_)*pow(abscissa, sizeExp_);
}

// ************************************************************************* //
