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

#include "univariateMomentSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMomentSet::univariateMomentSet
(
    const label nMoments,
    const supportType& support,
    const scalar smallM0,
    const scalar smallZeta,
    const scalar initValue,
    const label nAdditionalQuadraturePoints
)
:
    momentSet
    (
        nMoments,
        1,
        makeUnivariateMomentOrders(nMoments),
        List<supportType>(1, support),
        smallM0,
        smallZeta,
        initValue
    ),
    nAdditionalQuadraturePoints_(nAdditionalQuadraturePoints),
    alpha_(),
    beta_(),
    zeta_(),
    canonicalMoments_(),
    zRecurrence_(),
    status_()
{
    if (nAdditionalQuadraturePoints < 0)
    {
        FatalErrorInFunction
            << "The number of additional quadrature points must be positive."
            << nl
            << "    Number of additional quadrature points: "
            << nAdditionalQuadraturePoints << nl
            << exit(FatalError);
    }

    alpha_.setSize(nAlphaRecurrence(nAdditionalQuadraturePoints_), 0);
    beta_.setSize(nBetaRecurrence(nAdditionalQuadraturePoints_), 0);
    zeta_.setSize(nZetaRecurrence(nAdditionalQuadraturePoints_), 0);

    if (support == supportType::ZeroOne)
    {
        canonicalMoments_.setSize
        (
            nZetaRecurrence(nAdditionalQuadraturePoints_), 0
        );
    }
}


Foam::univariateMomentSet::univariateMomentSet
(
    const scalarList& m,
    const supportType& support,
    const scalar smallM0,
    const scalar smallZeta,
    const label nAdditionalQuadraturePoints
)
:
    momentSet
    (
        m,
        1,
        makeUnivariateMomentOrders(m.size()),
        List<supportType>(1, support),
        smallM0,
        smallZeta
    ),
    nAdditionalQuadraturePoints_(nAdditionalQuadraturePoints),
    alpha_(),
    beta_(),
    zeta_(),
    canonicalMoments_(),
    zRecurrence_(),
    status_()
{
    if (nAdditionalQuadraturePoints < 0)
    {
        FatalErrorInFunction
            << "The number of additional quadrature points must be positive."
            << nl
            << "    Number of additional quadrature points: "
            << nAdditionalQuadraturePoints << nl
            << exit(FatalError);
    }

    alpha_.setSize(nAlphaRecurrence(nAdditionalQuadraturePoints_), 0);
    beta_.setSize(nBetaRecurrence(nAdditionalQuadraturePoints_), 0);
    zeta_.setSize(nZetaRecurrence(nAdditionalQuadraturePoints_), 0);

    if (support == supportType::ZeroOne)
    {
        canonicalMoments_.setSize
        (
            nZetaRecurrence(nAdditionalQuadraturePoints_), 0
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::univariateMomentSet::nAlphaRecurrence
(
    const label nAdditionalQuadratureNodes
) const
{
    // With nMoments = 2p moments the recurrence relationship determines
    // alpha_0, ..., alpha_{p-1}; with nMoments = 2p + 1 moments it determines
    // the same coefficients, because alpha_p would require the moment of
    // order 2p + 1. In both cases the list holds nMoments/2 coefficients.
    return nMoments()/2 + nAdditionalQuadratureNodes;
}


Foam::label Foam::univariateMomentSet::nBetaRecurrence
(
    const label nAdditionalQuadratureNodes
) const
{
    return nMoments()/2 + 1 + nAdditionalQuadratureNodes;
}


Foam::label Foam::univariateMomentSet::nZetaRecurrence
(
    const label nAdditionalQuadratureNodes
) const
{
    // The moments determine nMoments - 1 values of zeta_k. The additional
    // quadrature nodes extend the chain to 2*nNodes - 1 values, which is what
    // GQMOM fills with a model number density function.
    return max
    (
        nMoments() - 1,
        2*(nMoments()/2 + nAdditionalQuadratureNodes) - 1
    );
}


Foam::labelListList Foam::univariateMomentSet::makeUnivariateMomentOrders
(
    const label nMoments
)
{
    labelListList mOrders(nMoments);

    for (label mI = 0; mI < nMoments; mI++)
    {
        mOrders[mI] = labelList(1, mI);
    }

    return mOrders;
}


void Foam::univariateMomentSet::resetStatus()
{
    status_ = realizabilityStatus();

    alpha_ = 0.0;
    beta_ = 0.0;
    zeta_ = 0.0;
    canonicalMoments_ = 0.0;
}


void Foam::univariateMomentSet::zetasToRecurrence
(
    const label nZeta,
    scalarList& alpha,
    scalarList& beta
) const
{
    if (supports_[0] == supportType::R)
    {
        FatalErrorInFunction
            << "The recurrence relationship cannot be recovered from the "
            << "zeta_k for measures with support over R." << nl
            << exit(FatalError);
    }

    if (nZeta < 1 || nZeta > zeta_.size())
    {
        FatalErrorInFunction
            << "The number of zeta_k is inconsistent with the size of the "
            << "zeta list." << nl
            << "    Number of zeta_k: " << nZeta << nl
            << "    Size of the zeta list: " << zeta_.size() << nl
            << exit(FatalError);
    }

    // Number of coefficients determined by nZeta values of the zeta chain.
    // alpha_i requires zeta_{2i} while beta_i only requires zeta_{2i-1}, so
    // with an even number of values one more beta is determined than alpha.
    const label nAlpha = (nZeta + 1)/2;
    const label nBeta = nZeta/2 + 1;

    if (nAlpha > alpha.size() || nBeta > beta.size())
    {
        FatalErrorInFunction
            << "The lists of the recurrence relationship are too small to "
            << "store the coefficients determined by the zeta_k." << nl
            << "    Number of zeta_k: " << nZeta << nl
            << "    Number of alpha coefficients: " << nAlpha
            << ", size of the alpha list: " << alpha.size() << nl
            << "    Number of beta coefficients: " << nBeta
            << ", size of the beta list: " << beta.size() << nl
            << exit(FatalError);
    }

    alpha[0] = zeta_[0];
    beta[0] = 1.0;

    for (label i = 1; i < nAlpha; i++)
    {
        alpha[i] = zeta_[2*i] + zeta_[2*i - 1];
    }

    for (label i = 1; i < nBeta; i++)
    {
        beta[i] = zeta_[2*i - 1]*zeta_[2*i - 2];
    }
}


void Foam::univariateMomentSet::update
(
    const scalarList& weights,
    const scalarList& abscissae
)
{
    // updateIntegerMoments resets the realizability status, because the
    // moments it recomputes are not necessarily realizable
    updateIntegerMoments(weights, abscissae);
}


void Foam::univariateMomentSet::updateIntegerMoments
(
    const scalarList& weights,
    const scalarList& abscissae
)
{
    // Recomputing all the moments (even if they originally were not realizable)
    // from quadrature (projection step).
    for (label momenti = 0; momenti < nMoments(); momenti++)
    {
        scalar& m = scalarList::operator[](momenti);

        m = Zero;

        forAll(weights, nodei)
        {
            m += weights[nodei]*pow(abscissae[nodei], momenti);
        }
    }

    resetStatus();
}


void Foam::univariateMomentSet::setSize(const label newSize)
{
    // Check that the new size is valid
    if (newSize < 2)
    {
        FatalErrorInFunction
            << "The new size of the moment set must be at least 2." << nl
            << "    New size: " << newSize << nl
            << exit(FatalError);
    }

    // Do not resize if the size is unchanged
    if (newSize == nMoments())
    {
        return;
    }

    labelListList newMomentOrders(makeUnivariateMomentOrders(newSize));

    // Resize the base moment set
    Foam::momentSet::setSize(newSize, newMomentOrders);

    // Resize the recurrence relationship and the zeta chain. These depend on
    // the new number of moments, so they are recomputed here.
    alpha_.setSize(nAlphaRecurrence(nAdditionalQuadraturePoints_), 0);
    beta_.setSize(nBetaRecurrence(nAdditionalQuadraturePoints_), 0);
    zeta_.setSize(nZetaRecurrence(nAdditionalQuadraturePoints_), 0);

    if (supports_[0] == supportType::ZeroOne)
    {
        canonicalMoments_.setSize
        (
            nZetaRecurrence(nAdditionalQuadraturePoints_), 0
        );
    }

    // Reset realizability status
    resetStatus();
}


void Foam::univariateMomentSet::resize(const label newSize)
{
    setSize(newSize);
}


// ************************************************************************* //
