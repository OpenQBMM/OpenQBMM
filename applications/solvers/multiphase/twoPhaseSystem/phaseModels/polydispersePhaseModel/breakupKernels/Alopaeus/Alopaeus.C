/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

#include "Alopaeus.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "phaseModel.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(Alopaeus, 0);

    addToRunTimeSelectionTable
    (
        breakupKernel,
        Alopaeus,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::Alopaeus::Alopaeus
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    breakupKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    C1_
    (
        dict.lookupOrDefault
        (
            "C1",
            dimensionedScalar("C1", dimless, 0.04)
        )
    ),
    C2_
    (
        dict.lookupOrDefault
        (
            "C2",
            dimensionedScalar("C2", dimless, 0.01)
        )
    ),
    epsilonf_
    (
        IOobject
        (
            "LuoSvendsen:epsilonf",
            fluid_.mesh().time().timeName(),
            fluid_.mesh()
        ),
        fluid_.mesh(),
        dimensionedScalar("zero", sqr(dimVelocity)/dimTime, 0.0)
    ),
    rhof_(fluid_.phase2().thermo().rho()()),
    muf_(fluid_.phase2().thermo().mu()()),
    rhob_(fluid_.phase1().thermo().rho()()),
    sigma_(fluid_.sigma())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernels::Alopaeus::~Alopaeus()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::populationBalanceSubModels::breakupKernels::Alopaeus::preUpdate()
{
    epsilonf_ = fluid_.phase2().turbulence().epsilon();
    epsilonf_.max(SMALL);
}


Foam::scalar
Foam::populationBalanceSubModels::breakupKernels::Alopaeus::Kb
(
    const scalar& d,
    const label celli,
    const label environment
) const
{
    return
        Cb_.value()*cbrt(epsilonf_[celli])
       *Foam::erfc
        (
            Foam::sqrt
            (
                C1_.value()*sigma_.value()
               /(
                    rhof_[celli]
                   *pow(epsilonf_[celli], 2.0/3.0)
                   *pow(d, 5.0/3.0)
                )
              + C2_.value()*muf_[celli]
               /(
                    sqrt(rhof_[celli]*rhob_[celli])
                   *cbrt(epsilonf_[celli])
                   *pow(d, 4.0/3.0)
                )
            )
        );
}

// ************************************************************************* //
