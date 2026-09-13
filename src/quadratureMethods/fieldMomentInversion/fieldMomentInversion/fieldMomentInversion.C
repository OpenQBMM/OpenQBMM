/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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

#include "fieldMomentInversion.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldMomentInversion, 0);
    defineRunTimeSelectionTable(fieldMomentInversion, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldMomentInversion::fieldMomentInversion
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
:
    mesh_(mesh),
    momentOrders_(momentOrders),
    nodeIndexes_(nodeIndexes),
    warnedSmallM0_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldMomentInversion::~fieldMomentInversion()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldMomentInversion::checkSmallM0(const volScalarField& m0) const
{
    if (warnedSmallM0_)
    {
        return;
    }

    // The zero-order moment carries the units of the weights, so smallM0 is
    // an absolute number that means different things in different cases.
    // With number densities it cuts nothing a real case holds; with the
    // volume fractions of small particles it can empty cells full of them,
    // and does so without saying. What can be told is how close the cut is
    // to the population: within a factor of this of the largest zero-order
    // moment in the domain, cells holding a millionth of that are taken as
    // empty, which is worth a warning.
    const scalar closeness = 1.0e6;

    // The boundaries are included: a case that injects through one starts
    // from an empty domain, with its population on the boundary until the
    // first steps carry it in.
    const scalar largestM0 = max(mag(m0)).value();

    if (largestM0 <= 0 || largestM0 >= closeness*smallM0())
    {
        return;
    }

    warnedSmallM0_ = true;

    WarningInFunction
        << "The largest zero-order moment of " << m0.name() << " is "
        << largestM0 << ", within a factor of " << closeness
        << " of smallM0, " << smallM0() << "." << nl
        << "    A cell whose zero-order moment is below smallM0 is taken as "
        << "empty, so part of this population may be discarded." << nl
        << "    If the weights are the volume fractions of small particles, "
        << "set smallM0 in the dictionary of the moment inversion below the "
        << "zero-order moment of the least populated cell that should be "
        << "kept." << nl
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
