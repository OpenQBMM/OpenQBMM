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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gaussMomentInversion::gaussMomentInversion
(
    univariateMomentSet& moments
)
:
    univariateMomentInversion(moments)
{
    if (moments_.nMoments() % 2 != 0)
    {
        FatalErrorInFunction
            << "The moment has an odd number of elements." << nl
            << "    Moment set: " << moments_
            << abort(FatalError);
    }

    calcNQuadratureNodes();

    weights_.setSize(nNodes_);
    abscissae_.setSize(nNodes_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gaussMomentInversion::~gaussMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gaussMomentInversion::calcNQuadratureNodes()
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
            << "The moment has size less or equal to 1." << nl
            << "    Moment set: " << moments_
            << abort(FatalError);
    }

    nNodes_ = nRealizableMoments/2;
}

// ************************************************************************* //
