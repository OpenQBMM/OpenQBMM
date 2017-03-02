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
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    scalarList& aRecurrence(moments.alphaRecurrence());
    scalarList& bRecurrence(moments.betaRecurrence());

    if (forceRadau_)
    {
        scalar p = minKnownAbscissa - aRecurrence[0];
        scalar pMinus1 = 1.0;
        scalar p1 = p;

        for (label i = 1; i < nNodes_ - 1; i++)
        {
            p = minKnownAbscissa - aRecurrence[0]*p1
                - moments.betaRecurrence()[i]*pMinus1;

            pMinus1 = p1;
            p1 = p;
        }

        aRecurrence[nNodes_ - 1] =
                minKnownAbscissa
              - moments.betaRecurrence()[nNodes_ - 1]*pMinus1/p;
    }
    else
    {
        scalar pLeft = minKnownAbscissa - aRecurrence[0];
        scalar pRight = maxKnownAbscissa - aRecurrence[0];

        scalar pMinus1Left = 1.0;
        scalar pMinus1Right = 1.0;

        scalar p1Left = pLeft;
        scalar p1Right = pRight;

        for (label i = 1; i < nNodes_ - 1; i++)
        {
            pLeft = (minKnownAbscissa - aRecurrence[i])*p1Left
                    - bRecurrence[i]*pMinus1Left;

            pRight = (maxKnownAbscissa - aRecurrence[i])*p1Right
                    - bRecurrence[i]*pMinus1Right;

            pMinus1Left = p1Left;
            pMinus1Right = p1Right;
            p1Left = pLeft;
            p1Right = pRight;
        }

        scalar d = pLeft*pMinus1Right - pRight*pMinus1Left;

        aRecurrence[nNodes_ - 1] =
                (minKnownAbscissa*pLeft*pMinus1Right
                - maxKnownAbscissa*pRight*pMinus1Left)/d;

        bRecurrence[nNodes_ - 1] =
                (maxKnownAbscissa - minKnownAbscissa)*pLeft*pRight/d;
    }
}

void Foam::gaussLobattoMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments,
    scalarList& weights,
    scalarList& abscissae
)
{
    if (moments.isDegenerate())
    {
        nNodes_ = 1;
        weights.setSize(nNodes_);
        abscissae.setSize(nNodes_);
        weights[0] = moments[0];
        abscissae[0] = 0.0;

        return;
    }

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

    weights.setSize(nNodes_);
    abscissae.setSize(nNodes_);
}

// ************************************************************************* //
