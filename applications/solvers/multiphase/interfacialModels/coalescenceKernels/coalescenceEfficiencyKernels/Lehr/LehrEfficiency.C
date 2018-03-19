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

#include "LehrEfficiency.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coalescenceEfficiencyKernels
{
    defineTypeNameAndDebug(Lehr, 0);
    addToRunTimeSelectionTable
    (
        coalescenceEfficiencyKernel,
        Lehr,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyKernels::Lehr::
Lehr
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalescenceEfficiencyKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    turbulence_(dict.lookupOrDefault("turbulence", true)),
    uCrit_
    (
        dimensionedScalar::lookupOrDefault
        (
            "uCrit",
            dict,
            dimVelocity,
            0.08
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceEfficiencyKernels::Lehr::
~Lehr()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::coalescenceEfficiencyKernels::Lehr::Pc
(
    const label nodei,
    const label nodej
) const
{
    const volScalarField& d1 = fluid_.phase1().ds(nodei);
    const volScalarField& d2 = fluid_.phase1().ds(nodej);

    if (turbulence_)
    {
        return max
        (
            uCrit_
           /max
            (
                1.414*cbrt(fluid_.phase2().turbulence().epsilon())
               *sqrt(pow(d1, 2/3) + pow(d2, 2/3)),
                mag(fluid_.phase1().Us(nodei) - fluid_.phase1().Us(nodej))
            ),
            1.0
        );
    }
    else
    {
        return max
        (
            uCrit_
           /max
            (
                mag(fluid_.phase1().Us(nodei) - fluid_.phase1().Us(nodej)),
                dimensionedScalar("smallVelocity", dimVelocity, 1e-10)
            ),
            1.0
        );
    }
}

// ************************************************************************* //
