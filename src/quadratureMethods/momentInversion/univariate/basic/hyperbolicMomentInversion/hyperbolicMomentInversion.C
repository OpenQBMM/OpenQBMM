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
    Copyright (C) 2019-2021 Alberto Passalacqua
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

#include "error.H"

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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hyperbolicMomentInversion::hyperbolicMomentInversion
(
    const dictionary& dict
)
:
    univariateMomentInversion(dict),
    etaMin_(dict.lookupOrDefault<scalar>("etaMin", 1.0e-10)),
    qMax_(dict.lookupOrDefault<scalar>("qMax", 30.0)),
    smallNegRealizability_
    (
        dict.lookupOrDefault<scalar>
        (
            "smallNegRealizability",
            -1.0e-6
        )
    )
{
    nInvertibleMoments_ = 5;
    nNodes_ = 3;
    weights_.setSize(nNodes_, Zero);
    abscissae_.setSize(nNodes_, Zero);
}


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
    NotImplemented;
    return;
}

void Foam::hyperbolicMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    nInvertibleMoments_ = 5;
    nNodes_ = 3;
}

void Foam::hyperbolicMomentInversion::invert
(
    univariateMomentSet& moments,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    if (moments[0] < SMALL)
    {
        weights_[0] = 0.0;
        weights_[1] = moments[0];
        weights_[2] = 0.0;

        abscissae_[0] = 0.0;
        abscissae_[1] = 0.0;
        abscissae_[2] = 0.0;

        return;
    }

    // Compute normalized moments
    scalarList normalizedMoments(moments);

    forAll(normalizedMoments, mi)
    {
        normalizedMoments[mi] /= moments[0];
    }

    scalar meanVelocity = normalizedMoments[1];

    // Compute central moments
    scalarList centralMoments(normalizedMoments);

    centralMoments[0] = 1.0;
    centralMoments[2] -= sqr(meanVelocity);

    centralMoments[3] -=
        3.0*meanVelocity*normalizedMoments[2]
      - 2.0*pow3(meanVelocity);

    centralMoments[4] -=
        4.0*meanVelocity*normalizedMoments[3]
      - 6.0*sqr(meanVelocity)*normalizedMoments[2]
      + 3.0*pow4(meanVelocity);

    // Compute realizability condition
    scalar realizability =
        centralMoments[2]*centralMoments[4]
      - pow3(centralMoments[2])
      - sqr(centralMoments[3]);

    // Manage unrealizable cases
    if (centralMoments[2] < 0.0)
    {
        if (centralMoments[2] < -1e-10)
        {
            WarningInFunction
                << "Second-order central moment is negative. C2 = "
                << centralMoments[2]
                << endl;
        }

        for (label ci = 2; ci < nInvertibleMoments_; ci++)
        {
            centralMoments[ci] = 0.0;
        }
    }
    else if (realizability < 0)
    {
        if (centralMoments[2] >= etaMin_)
        {
            scalar c2 = centralMoments[2];
            scalar sqrC2 = sqr(c2);
            scalar sqrtC2 = sqrt(c2);

            scalar q = centralMoments[3]/(c2*sqrtC2);
            scalar eta = centralMoments[4]/sqrC2;

            if (mag(q) > SMALL)
            {
                scalar slope = (eta - 3.0)/q;
                scalar sqrtDet = sqrt(8.0 + sqr(slope));

                if (q > 0.0)
                {
                    q = (slope + sqrtDet)/2.0;
                }
                else
                {
                    q = (slope - sqrtDet)/2.0;
                }
            }
            else
            {
                q = 0.0;
            }

            eta = 1.0 + sqr(q);
            centralMoments[3] = q*c2*sqrtC2;
            centralMoments[4] = eta*sqrC2;

            if (realizability < smallNegRealizability_ && debug)
            {
                WarningInFunction
                    << "Fourth-order central moment is too SMALL."
                    << " Realizability = " << realizability << nl
                    << endl;
            }
        }
        else
        {
            centralMoments[3] = 0.0;
            centralMoments[4] = sqr(centralMoments[2]);
        }
    }

    // Find scaled parameters
    scalar c2 = centralMoments[2];
    scalar sqrC2 = sqr(c2);
    scalar sqrtC2 = sqrt(c2);
    scalar q = 0.0;
    scalar eta = 0.0;

    if (c2 >= etaMin_)
    {
        q = centralMoments[3]/(c2*sqrtC2);
        eta = centralMoments[4]/sqrC2;
    }
    else
    {
        q = 0.0;
        eta = 1.0;
    }

    // Bind skewness
    if (sqr(q) > sqr(qMax_))
    {
        scalar slope = (eta - 3.0)/q;
        q = qMax_*sign(q);
        eta = 3.0 + q*slope;
        realizability = eta - sqr(q) - 1.0;

        if (realizability < 0.0)
        {
            eta = 1.0 + sqr(q);
        }
    }

    scalar sqrQ = sqr(q);
    scalar sqrtZ = sqrt(4.0*eta - 3.0*sqrQ);
    scalar invSqrtZ = 1.0/sqrtZ;

    // Compute unscaled abscissae
    abscissae_[0] = 0.5*(q - sqrtZ);
    abscissae_[1] = 0.0;
    abscissae_[2] = 0.5*(q + sqrtZ);

    scalar abscissaeProd = max(- abscissae_[0]*abscissae_[2], 1.0 + SMALL);

    // Compute weights
    weights_[0] = - invSqrtZ/abscissae_[0];
    weights_[1] = 1.0 - 1.0/abscissaeProd;
    weights_[2] = invSqrtZ/abscissae_[2];

    // Sum weights
    scalar weightSum = 0.0;
    forAll(weights_, wi)
    {
        weightSum += weights_[wi];
    }

    // Rescale weights
    forAll(weights_, wi)
    {
        weights_[wi] /= weightSum;
    }

    scalar scaleFactor = 0.0;

    forAll(weights_, wi)
    {
        scaleFactor += weights_[wi]*sqr(abscissae_[wi])/weightSum;
    }

    // Rescale abscissae
    forAll(abscissae_, ai)
    {
        abscissae_[ai] *= sqrtC2/scaleFactor;
    }

    forAll(weights_, wi)
    {
        weights_[wi] *= moments[0];

        if (weights_[wi] < 0.0)
        {
            FatalErrorInFunction
                << "Negative weight in hyperbolic moment inversion."
                << exit(FatalError);
        }

        abscissae_[wi] += meanVelocity;
    }
}

// ************************************************************************* //
