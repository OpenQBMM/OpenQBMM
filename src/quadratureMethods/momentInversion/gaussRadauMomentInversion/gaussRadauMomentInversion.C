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

#include "gaussRadauMomentInversion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gaussRadauMomentInversion::gaussRadauMomentInversion
(
    univariateMomentSet& moments,
    scalar knownAbscissa
)
:
    univariateMomentInversion(moments),
    forceGauss_(false),
    knownAbscissa_(knownAbscissa)
{
    if (moments_.nMoments() % 2 == 0)
    {
        FatalErrorInFunction
            << "The moment has an even number of elements." << nl
            << "    Moment set: " << moments_
            << abort(FatalError);
    }

    calcNQuadratureNodes();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gaussRadauMomentInversion::~gaussRadauMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gaussRadauMomentInversion::updateRecurrenceRelation()
{
    scalarList& aRecurrence(moments_.alphaRecurrence());

    if (!forceGauss_)
    {
        scalar p = knownAbscissa_ - aRecurrence[0];
        scalar pMinus1 = 1.0;
        scalar p1 = p;

        for (label i = 1; i < nNodes_ - 1; i++)
        {
            p = knownAbscissa_ - aRecurrence[0]*p1
                - moments_.betaRecurrence()[i]*pMinus1;

            pMinus1 = p1;
            p1 = p;
        }

        aRecurrence[nNodes_ - 1] =
            knownAbscissa_
          - moments_.betaRecurrence()[nNodes_ - 1]*pMinus1/p;
    }
}

void Foam::gaussRadauMomentInversion::calcNQuadratureNodes()
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
            forceGauss_ = true;
            nNodes_ = nInvertibleMoments_/2;
        }
        else
        {
            nInvertibleMoments_ = nRealizableMoments;
            forceGauss_ = false;
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

void Foam::gaussRadauMomentInversion::invert()
{
    if (!inverted_)
    {
        calcNQuadratureNodes();
        updateRecurrenceRelation();

        univariateMomentInversion::invert();
    }
}


// ************************************************************************* //
