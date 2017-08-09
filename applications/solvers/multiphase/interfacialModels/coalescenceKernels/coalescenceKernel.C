/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
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

#include "coalescenceKernel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalescenceKernel::coalescenceKernel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    Ca_
    (
        dict.lookupOrDefault
        (
            "Ca",
            dimensionedScalar("one", dimless, 1.0)
        )
    ),
    frequency_(coalescenceFrequencyKernel::New(dict, mesh)),
    efficiency_(coalescenceEfficiencyKernel::New(dict, mesh))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coalescenceKernel::~coalescenceKernel()
{}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::coalescenceKernel::Ka
(
    const label nodei,
    const label nodej
) const
{
    return
        Ca_
       *frequency_->omega(nodei, nodej)
       *efficiency_->Pc(nodei, nodej);
}


// ************************************************************************* //
