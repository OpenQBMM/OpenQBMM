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
    Copyright (C) 2019 Alberto Passalacqua
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

#include "gaussRadauMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gaussRadauMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentInversion,
        gaussRadauMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gaussRadauMomentInversion::gaussRadauMomentInversion
(
    const dictionary& dict
)
:
    univariateMomentInversion(dict),
    forceGauss_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gaussRadauMomentInversion::~gaussRadauMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gaussRadauMomentInversion::correctRecurrence
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    if (!forceGauss_)
    {
        scalar p = minKnownAbscissa - alpha[0];
        scalar pMinus1 = 1.0;
        scalar p1 = p;

        for (label i = 1; i < nNodes_ - 1; i++)
        {
            p = (minKnownAbscissa - alpha[i])*p1 - beta[i]*pMinus1;

            pMinus1 = p1;
            p1 = p;
        }

        alpha[nNodes_ - 1] =
            minKnownAbscissa - beta[nNodes_ - 1]*pMinus1/p;
    }
}

void Foam::gaussRadauMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    label nRealizableMoments = moments.nRealizableMoments();

    if (nRealizableMoments >= 2)
    {
        if (nRealizableMoments % 2 == 0)
        {
            nInvertibleMoments_ = nRealizableMoments;
            forceGauss_ = true;
            nNodes_ = nInvertibleMoments_/2;
        }
        else
        {
            nInvertibleMoments_ = nRealizableMoments;
            forceGauss_ = false;
            nNodes_ = (nInvertibleMoments_ - 1)/2 + 1;
        }
    }
    else
    {
        FatalErrorInFunction
            << "The moment has size less or equal to 1." << nl
            << "    Moment set: " << moments
            << abort(FatalError);
    }

    abscissae_.setSize(nNodes_);
    weights_.setSize(nNodes_);
}

// ************************************************************************* //
