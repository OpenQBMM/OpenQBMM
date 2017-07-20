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

#include "Chesters.H"
#include "addToRunTimeSelectionTable.H"
#include "fundamentalConstants.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace coalesenceEffeciencyKernels
{
    defineTypeNameAndDebug(Chesters, 0);

    addToRunTimeSelectionTable
    (
        coalesenceEffeciencyKernel,
        Chesters,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::coalesenceEffeciencyKernels::
Chesters::Chesters
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coalesenceEffeciencyKernel(dict, mesh),
    fluid_(mesh.lookupObject<twoPhaseSystem>("phaseProperties")),
    Ceff_(dict.lookup("Ceff"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::coalesenceEffeciencyKernels::
Chesters::~Chesters()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalanceSubModels::coalesenceEffeciencyKernels::
Chesters::Pc
(
    const volScalarField& d1,
    const volScalarField& d2
) const
{
    const phasePair& pair = fluid_.pair();
    volScalarField We = pair.We();
    const dimensionedScalar& sigma = fluid_.sigma();

    return
        Foam::exp
        (
          - Ceff_
           *sqrt
            (
                nu*epsilon*sqr(rho/sigma)
               *pow4
                (
                    d1*d2
                   /max(d1 + d2, dimensionedScalar("zero", dimLength, SMALL))
                )
            )
        );
}

// ************************************************************************* //
