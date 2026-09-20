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
    Copyright (C) 2019-2026 Alberto Passalacqua
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

#include "univariateMomentInversion.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(univariateMomentInversion, 0);
    defineRunTimeSelectionTable(univariateMomentInversion, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMomentInversion::univariateMomentInversion
(
    const dictionary& dict,
    const label nMaxNodes
)
:
    smallM0_(dict.lookupOrDefault<scalar>("smallM0", SMALL)),
    smallZeta_(dict.lookupOrDefault<scalar>("smallZeta", 0.0)),
    boundaryTolerance_
    (
        dict.lookupOrDefault<scalar>("boundaryTolerance", 1e-12)
    ),
    nInvertibleMoments_(),
    nNodes_(nMaxNodes),
    abscissae_(),
    weights_(),
    golubWelsch_(),
    alpha_(),
    beta_()
{
    if (smallZeta_ < 0.0)
    {
        FatalErrorInFunction
            << "The value of smallZeta must be positive or null."
            << exit(FatalError);
    }

    if (boundaryTolerance_ < 0 || boundaryTolerance_ >= 0.5)
    {
        FatalIOErrorInFunction(dict)
            << "The value of boundaryTolerance must be in [0, 0.5)." << nl
            << "    boundaryTolerance: " << boundaryTolerance_ << nl
            << exit(FatalIOError);
    }

    if (smallZeta_ > 0)
    {
        WarningInFunction
            << "The value of smallZeta is larger than zero. " << endl
            << "This may lead to the exclusion of valid moment vectors." << endl
            << endl
            << "smallZeta = " << smallZeta_
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::univariateMomentInversion::nAdditionalQuadraturePoints
(
    const label nMoments
) const
{
    // A quadrature of nMoments/2 nodes is what the moment set is sized for,
    // so a quadrature that does not go beyond it needs no additional room
    return 0;
}


Foam::scalar Foam::univariateMomentInversion::orthogonalPolynomial
(
    const scalarList& alpha,
    const scalarList& beta,
    const scalar x,
    scalar& pMinus1
) const
{
    // Three-term recurrence of the monic orthogonal polynomials:
    // p_{k+1}(x) = (x - alpha_k) p_k(x) - beta_k p_{k-1}(x), with p_0 = 1
    scalar p = x - alpha[0];

    pMinus1 = 1.0;

    for (label i = 1; i < nNodes_ - 1; i++)
    {
        const scalar pNext = (x - alpha[i])*p - beta[i]*pMinus1;

        pMinus1 = p;
        p = pNext;
    }

    return p;
}


void Foam::univariateMomentInversion::correctedRecurrence
(
    univariateMomentSet& moments
)
{
    const scalarList& momentAlpha = moments.alphaRecurrence();
    const scalarList& momentBeta = moments.betaRecurrence();

    if (nNodes_ > momentAlpha.size() || nNodes_ > momentBeta.size())
    {
        FatalErrorInFunction
            << "The recurrence relationship of the moment set is too short "
            << "for the requested number of quadrature nodes." << nl
            << "    Number of quadrature nodes: " << nNodes_ << nl
            << "    Size of the alpha list: " << momentAlpha.size() << nl
            << "    Size of the beta list: " << momentBeta.size() << nl
            << "    Moment set: " << moments << nl
            << exit(FatalError);
    }

    // correctRecurrence builds the coefficients of the quadrature, leaving
    // the ones of the moment set untouched, so they are copied into the work
    // lists. Both are members, so the copy does not allocate once the
    // inversion has run on the first cell.
    alpha_ = momentAlpha;
    beta_ = momentBeta;

    correctRecurrence
    (
        moments,
        alpha_,
        beta_
    );
}

void Foam::univariateMomentInversion::invert
(
    univariateMomentSet& moments
)
{
    // The moments are read through the const accessor, so that the
    // realizability check performed by the caller is not invalidated and
    // recomputed by calcNQuadratureNodes below
    const scalar m0 = moments.moment(0);

    if (moments.isDegenerate())
    {
        nNodes_ = 1;
        weights_.setSize(nNodes_);
        abscissae_.setSize(nNodes_);
        weights_[0] = m0;
        abscissae_[0] = 0.0;

        return;
    }

    if (m0 < smallM0_)
    {
        nNodes_ = 0;

        weights_.clear();
        abscissae_.clear();

        return;
    }

    if (invertOnBoundary(moments))
    {
        return;
    }

    calcNQuadratureNodes(moments);

    if (nInvertibleMoments_ == 2)
    {
        weights_[0] = m0;
        abscissae_[0] = moments.moment(1)/m0;

        return;
    }

    correctedRecurrence(moments);
    calcQuadrature(moments);
}

bool Foam::univariateMomentInversion::invertOnBoundary
(
    univariateMomentSet& moments
)
{
    if (moments.support() != supportType::ZeroOne)
    {
        return false;
    }

    // The tolerance is the one of this inversion, which the moment set,
    // built by whoever owns it, does not know
    moments.boundaryTolerance(boundaryTolerance_);

    const label N = moments.boundaryOrder(false);

    if (N < 2)
    {
        // In the interior of the moment space, or all of the measure at
        // one, which the mean gives
        return false;
    }

    // The zeta_k of the measure, zeta[k] being zeta_k with zeta_0 = 0, up
    // to the first that is zero
    const scalarList& setZeta = moments.zetas();
    const bool atOne = moments.boundaryAtOne(false);
    const label nZeta = atOne ? N : N - 1;

    scalarList zeta(nZeta + 2, Zero);

    for (label k = 1; k <= nZeta; k++)
    {
        zeta[k] = setZeta[k - 1];
    }

    if (atOne)
    {
        // zeta_N = p_N (1 - p_(N-1)) with p_N at one
        zeta[N] = 1.0 - moments.canonicalMoments()[N - 2];
    }

    // The first zeta that is zero is zeta_M
    const label M = nZeta + 1;

    nNodes_ = (M + 1)/2;
    nInvertibleMoments_ = N;

    weights_.setSize(nNodes_);
    abscissae_.setSize(nNodes_);
    alpha_.setSize(nNodes_);
    beta_.setSize(nNodes_);

    auto zetaAt = [&zeta](const label k)
    {
        return k < zeta.size() ? zeta[k] : scalar(0);
    };

    for (label j = 0; j < nNodes_; j++)
    {
        alpha_[j] = zetaAt(2*j) + zetaAt(2*j + 1);
        beta_[j] = j > 0 ? zetaAt(2*j - 1)*zetaAt(2*j) : scalar(1);
    }

    calcQuadrature(moments);

    // The ends of the interval are nodes to round-off
    forAll(abscissae_, nodei)
    {
        abscissae_[nodei] = min(max(abscissae_[nodei], scalar(0)), scalar(1));
        weights_[nodei] = max(weights_[nodei], scalar(0));
    }

    return true;
}


void Foam::univariateMomentInversion::calcQuadrature
(
    const univariateMomentSet& moments
)
{
    golubWelsch_.quadrature
    (
        alpha_,
        beta_,
        nNodes_,
        moments.moment(0),
        abscissae_,
        weights_
    );
}

// ************************************************************************* //
