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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gaussLobattoMomentInversion::gaussLobattoMomentInversion
(
    univariateMomentSet& moments,
    scalar minKnownAbscissa,
    scalar maxKnownAbscissa
)
:
    univariateMomentInversion(moments),
    forceRadau_(false),
    minKnownAbscissa_(minKnownAbscissa),
    maxKnownAbscissa_(maxKnownAbscissa)
{
    if (minKnownAbscissa >= maxKnownAbscissa)
    {
        FatalErrorInFunction
            << "The interval of integration is not correctly specified." << nl
            << "    Min. abscissa: " << minKnownAbscissa << nl
            << "    Max. abscissa: " << maxKnownAbscissa
            << abort(FatalError);
    }

    if (moments_.nMoments() % 2 != 0)
    {
        FatalErrorInFunction
            << "The moment has an odd number of elements." << nl
            << "    Moment set: " << moments_
            << abort(FatalError);
    }

    calcNQuadratureNodes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gaussLobattoMomentInversion::~gaussLobattoMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gaussLobattoMomentInversion::updateRecurrenceRelation()
{
    scalarList& aRecurrence(moments_.alphaRecurrence());
    scalarList& bRecurrence(moments_.betaRecurrence());

    Info << "a " << aRecurrence[nNodes_ - 1] << endl;
    Info << "b " << bRecurrence[nNodes_ - 1] << endl;

    if (forceRadau_)
    {
        scalar p = minKnownAbscissa_ - aRecurrence[0];
        scalar pMinus1 = 1.0;
        scalar p1 = p;

        for (label i = 1; i < nNodes_ - 1; i++)
        {
            p = minKnownAbscissa_ - aRecurrence[0]*p1
                - moments_.betaRecurrence()[i]*pMinus1;

            pMinus1 = p1;
            p1 = p;
        }

        aRecurrence[nNodes_ - 1] =
                minKnownAbscissa_
              - moments_.betaRecurrence()[nNodes_ - 1]*pMinus1/p;
    }
    else
    {
        scalar pLeft = minKnownAbscissa_ - aRecurrence[0];
        scalar pRight = maxKnownAbscissa_ - aRecurrence[0];

        scalar pMinus1Left = 1.0;
        scalar pMinus1Right = 1.0;

        scalar p1Left = pLeft;
        scalar p1Right = pRight;

        for (label i = 1; i < nNodes_ - 1; i++)
        {
            pLeft = (minKnownAbscissa_ - aRecurrence[i])*p1Left
                    - bRecurrence[i]*pMinus1Left;

            pRight = (maxKnownAbscissa_ - aRecurrence[i])*p1Right
                    - bRecurrence[i]*pMinus1Right;

            pMinus1Left = p1Left;
            pMinus1Right = p1Right;
            p1Left = pLeft;
            p1Right = pRight;
        }

        scalar d = pLeft*pMinus1Right - pRight*pMinus1Left;

        Info << "det" << d << endl;

        aRecurrence[nNodes_ - 1] =
                (minKnownAbscissa_*pLeft*pMinus1Right
                - maxKnownAbscissa_*pRight*pMinus1Left)/d;

        bRecurrence[nNodes_ - 1] =
                (maxKnownAbscissa_ - minKnownAbscissa_)*pLeft*pRight/d;
    }

    Info << "a " << aRecurrence[nNodes_ - 1] << endl;
    Info << "b " << bRecurrence[nNodes_ - 1] << endl;
}


void Foam::gaussLobattoMomentInversion::calcNQuadratureNodes()
{
    if (moments_.isDegenerate())
    {
        nNodes_ = 1;
        weights_.setSize(nNodes_);
        abscissae_.setSize(nNodes_);
        weights_[0] = moments_[0];
        abscissae_[0] = 0.0;
        inverted_ = true;

        return;
    }

    label nRealizableMoments = moments_.nRealizableMoments();

    if (nRealizableMoments >= 2)
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
            << "The moment has size less or equal to 1." << nl
            << "    Moment set: " << moments_
            << abort(FatalError);
    }

    weights_.setSize(nNodes_);
    abscissae_.setSize(nNodes_);
}

void Foam::gaussLobattoMomentInversion::invert()
{
    if (!inverted_)
    {
        calcNQuadratureNodes();
        updateRecurrenceRelation();

        univariateMomentInversion::invert();
    }
}

// ************************************************************************* //
