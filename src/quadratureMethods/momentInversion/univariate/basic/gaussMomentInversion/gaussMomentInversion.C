/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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

#include "gaussMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gaussMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentInversion,
        gaussMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gaussMomentInversion::gaussMomentInversion
(
    const dictionary& dict
)
:
    univariateMomentInversion(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gaussMomentInversion::~gaussMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gaussMomentInversion::correctRecurrence
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    return; // No correction needed for Gauss quadrature
}

void Foam::gaussMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    label nRealizableMoments = moments.nRealizableMoments();

    if (nRealizableMoments >= 2)
    {
        if (nRealizableMoments % 2 != 0)
        {
            nInvertibleMoments_ = nRealizableMoments - 1;
        }
        else
        {
            nInvertibleMoments_ = nRealizableMoments;
        }
    }
    else
    {
        FatalErrorInFunction
            << "The moment set has size less or equal to 1." << nl
            << "    Moment set: " << moments
            << abort(FatalError);
    }

    nNodes_ = nInvertibleMoments_/2;

    abscissae_.setSize(nNodes_);
    weights_.setSize(nNodes_);
}

// ************************************************************************* //
