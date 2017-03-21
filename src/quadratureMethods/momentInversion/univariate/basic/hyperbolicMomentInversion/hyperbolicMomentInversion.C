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

#include "hyperbolicMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hyperbolicMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentInversion,
        hyperbolicMomentInversion,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::hyperbolicMomentInversion::calculateCentralMoments
(
    const scalarList& moments
)
{
    centralMoments_[0] = scalar(1);
    centralMoments_[1] = scalar(0);

    if (nNodes_ == 1)
    {
        return;
    }
    if (nNodes_ > 1)
    {
        centralMoments_[2] =
            (
                moments[0]*moments[2]
              - Foam::sqr(moments[1])
            )/Foam::sqr(moments[0]);

        return;
    }
    if (nNodes_ > 2)
    {
        centralMoments_[3] =
            (
                Foam::sqr(moments[0])*moments[3]
              - 3.0*moments[0]*moments[1]*moments[2]
              + 2.0*Foam::pow3(moments[1])
            )/Foam::pow3(moments[0]);

        centralMoments_[4] =
            (
                Foam::pow3(moments[0])*moments[4]
              - 4.0*Foam::sqr(moments[0])*moments[1]*moments[3]
              + 6.0*moments[0]*Foam::sqr(moments[1])*moments[2]
              - 3.0*Foam::pow4(moments[1])
            )/Foam::pow4(moments[0]);

        return;
    }
    if (nNodes_ > 3)
    {
        FatalErrorInFunction
            << "Hyperbolic moment inversion can handle a maximum of " << nl
            << "3 nodes. " << nNodes_ << " were selected."
            << abort(FatalError);
    }
}

void Foam::hyperbolicMomentInversion::unscaleAbscissae
(
    const scalarList& moments
)
{
    scalar m0 = moments[0];
    scalar uBar = moments[1]/moments[0];

    forAll(abscissae_, nodei)
    {
        abscissae_[nodei] += uBar;
        weights_[nodei] *= m0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hyperbolicMomentInversion::hyperbolicMomentInversion
(
    const dictionary& dict
)
:
    univariateMomentInversion(dict),
    centralMoments_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hyperbolicMomentInversion::~hyperbolicMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hyperbolicMomentInversion::correctRecurrence
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


void Foam::hyperbolicMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    label nRealizableMoments = moments.nRealizableMoments();
    if (nRealizableMoments >= 2)
    {
        if ((nRealizableMoments + 1) % 2 != 0)
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

    nNodes_ = (nInvertibleMoments_ + 1)/2;

    if (nNodes_ > 3)
    {
        FatalErrorInFunction
            << "Hyperbolic moment inversion can handle a maximum of " << nl
            << "3 nodes. " << nNodes_ << " were selected."
            << abort(FatalError);
    }

    abscissae_.setSize(nNodes_);
    weights_.setSize(nNodes_);
}


void Foam::hyperbolicMomentInversion::invert
(
    univariateMomentSet& moments,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    calcNQuadratureNodes(moments);

    centralMoments_ = scalarList(2*nNodes_, scalar(0));
    calculateCentralMoments(moments);

    if (nNodes_ == 1)
    {
        weights_[0] = moments[0];
        abscissae_[0] = moments[1]/moments[0];

        return;
    }
    if (nNodes_ == 2)
    {
        //  Set weights
        weights_[0] = 0.5;
        weights_[1] = 0.5;

        //  Set abscissae
        abscissae_[0] = Foam::sqrt(centralMoments_[2]);
        abscissae_[1] = -1*abscissae_[0];

        //  Return to normal moment components
        unscaleAbscissae(moments);

        return;
    }
    if (nNodes_ == 3)
    {
        //  Calculate intermediate variables
        scalar q =
            centralMoments_[3]
           /(centralMoments_[2]*Foam::sqrt(centralMoments_[2]));

        scalar eta =
            centralMoments_[4]/Foam::sqr(centralMoments_[2]);

        //  Set abscissae
        abscissae_[0] =
            0.5*Foam::sqrt(centralMoments_[2])
           *(q - Foam::sqrt(4.0*eta - 3.0*Foam::sqr(q)));

        abscissae_[1] = scalar(0);

        abscissae_[2] =
            0.5*Foam::sqrt(centralMoments_[2])
           *(q + Foam::sqrt(4.0*eta - 3.0*Foam::sqr(q)));

        //  Set weights
        weights_[0] =
           -Foam::sqr(centralMoments_[2])
           /(
                abscissae_[0]
               *Foam::sqrt
                (
                    4.0*centralMoments_[2]*centralMoments_[4]
                  - 3.0*Foam::sqr(centralMoments_[3])
                )
            );

        weights_[1] =
            1.0
          + centralMoments_[2]
           /(abscissae_[0]*abscissae_[2]);

        weights_[2] =
            Foam::sqr(centralMoments_[2])
           /(
                abscissae_[2]
               *Foam::sqrt
                (
                    4.0*centralMoments_[2]*centralMoments_[4]
                  - 3.0*Foam::sqr(centralMoments_[3])
                )
            );

        //  Return to normal moment components
        unscaleAbscissae(moments);

        return;
    }
    if (nNodes_ > 3)
    {
        FatalErrorInFunction
            << "Hyperbolic moment inversion can handle a maximum of " << nl
            << "3 nodes. " << nNodes_ << " were selected."
            << abort(FatalError);
    }
}
// ************************************************************************* //
