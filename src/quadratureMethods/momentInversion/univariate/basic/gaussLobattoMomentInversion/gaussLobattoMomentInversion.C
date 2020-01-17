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

#include "gaussLobattoMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gaussLobattoMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentInversion,
        gaussLobattoMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gaussLobattoMomentInversion::gaussLobattoMomentInversion
(
    const dictionary& dict
)
:
    univariateMomentInversion(dict),
    forceRadau_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gaussLobattoMomentInversion::~gaussLobattoMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gaussLobattoMomentInversion::correctRecurrence
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    if (forceRadau_)
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
    else
    {
        scalar pLeft = minKnownAbscissa - alpha[0];
        scalar pRight = maxKnownAbscissa - alpha[0];

        scalar pMinus1Left = 1.0;
        scalar pMinus1Right = 1.0;

        scalar p1Left = pLeft;
        scalar p1Right = pRight;

        for (label i = 1; i < nNodes_ - 1; i++)
        {
            pLeft = (minKnownAbscissa - alpha[i])*p1Left
                    - beta[i]*pMinus1Left;

            pRight = (maxKnownAbscissa - alpha[i])*p1Right
                    - beta[i]*pMinus1Right;

            pMinus1Left = p1Left;
            pMinus1Right = p1Right;
            p1Left = pLeft;
            p1Right = pRight;
        }

        scalar d = pLeft*pMinus1Right - pRight*pMinus1Left;

        alpha[nNodes_ - 1] =
                (minKnownAbscissa*pLeft*pMinus1Right
                - maxKnownAbscissa*pRight*pMinus1Left)/d;

        beta[nNodes_ - 1] =
                (maxKnownAbscissa - minKnownAbscissa)*pLeft*pRight/d;
    }
}

void Foam::gaussLobattoMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    label nRealizableMoments = moments.nRealizableMoments();

    if (nRealizableMoments > 2)
    {
        if (nRealizableMoments % 2 == 0)
        {
            nInvertibleMoments_ = nRealizableMoments;
            forceRadau_ = false;
            nNodes_ = nInvertibleMoments_/2 + 1;
        }
        else
        {
            nInvertibleMoments_ = nRealizableMoments;
            forceRadau_ = true;
            nNodes_ = nInvertibleMoments_/2 + 1;
        }
    }
    else
    {
        FatalErrorInFunction
            << "The moment has size less or equal to 2." << nl
            << "    Moment set: " << moments
            << abort(FatalError);
    }

    abscissae_.setSize(nNodes_);
    weights_.setSize(nNodes_);
}

// ************************************************************************* //
