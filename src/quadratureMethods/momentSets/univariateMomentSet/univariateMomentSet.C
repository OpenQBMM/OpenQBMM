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
    Copyright (C) 2019-2023 Alberto Passalacqua
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
    const word& support,
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
        support,
        smallM0,
        smallZeta,
        initValue
    ),
    alpha_(),
    beta_(),
    zeta_(nMoments_ - 1),
    canonicalMoments_(),
    negativeZeta_(0),
    degenerate_(false),
    fullyRealizable_(true),
    subsetRealizable_(true),
    onMomentSpaceBoundary_(false),
    nRealizableMoments_(0),
    realizabilityChecked_(false)
{
    if (support_ != "R" && support_ != "RPlus" && support_ != "01")
    {
        FatalErrorInFunction
            << "The specified support is invalid." << nl
            << "    Valid supports are: R, RPlus and 01."
            << abort(FatalError);
    }

    if (nAdditionalQuadraturePoints < 0)
    {
        FatalErrorInFunction
            << "The number of additional quadrature points must be positive."
            << abort(FatalError);
    }
   
    label nAlpha = nAlphaRecurrence(nAdditionalQuadraturePoints);
    label nBeta = nBetaRecurrence(nAdditionalQuadraturePoints);;

    alpha_.setSize(nAlpha, 0);
    beta_.setSize(nBeta, 0);

    if (support_ == "01")
    {
        canonicalMoments_.setSize(nMoments_ - 1, 0);
    }
}

Foam::univariateMomentSet::univariateMomentSet
(
    const scalarList& m,
    const word& support,
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
        support,
        smallM0,
        smallZeta
    ),
    alpha_(),
    beta_(),
    zeta_(nMoments_ - 1),
    negativeZeta_(0),
    degenerate_(false),
    fullyRealizable_(true),
    subsetRealizable_(true),
    onMomentSpaceBoundary_(false),
    nRealizableMoments_(0),
    realizabilityChecked_(false)
{
    if (support_ != "R" && support_ != "RPlus" && support_ != "01")
    {
        FatalErrorInFunction
            << "The specified support is invalid." << nl
            << "    Valid supports are: R, RPlus and 01."
            << abort(FatalError);
    }

    if (nAdditionalQuadraturePoints < 0)
    {
        FatalErrorInFunction
            << "The specified number of fixed points must be positive." << nl
            << abort(FatalError);
    }

    label nAlpha = nAlphaRecurrence(nAdditionalQuadraturePoints);
    label nBeta = nBetaRecurrence(nAdditionalQuadraturePoints);

    alpha_.setSize(nAlpha, 0);
    beta_.setSize(nBeta, 0);

    if (support_ == "01")
    {
        canonicalMoments_.setSize(nMoments_ - 1, 0);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateMomentSet::~univariateMomentSet()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::univariateMomentSet::nAlphaRecurrence
(
    const label& nAdditionalQuadraturePoints
)
{
    return label((nMoments_ - 2)/2) + 1 + nAdditionalQuadraturePoints;
}


Foam::label Foam::univariateMomentSet::nBetaRecurrence
(
    const label& nAdditionalQuadraturePoints
)
{
    return label(nMoments_/2) + 1 + nAdditionalQuadraturePoints;
}


void Foam::univariateMomentSet::checkCanonicalMoments
(
    const scalarList& zeta,
    const label nZeta
)
{
    canonicalMoments_ = 0.0;
    canonicalMoments_[0] = zeta[0];

    if (mag(canonicalMoments_[0] - 1.0) <= smallZeta_)
    {
        nRealizableMoments_ = 2;
        onMomentSpaceBoundary_ = true;

        return;
    }

    for (label zetai = 1; zetai < nZeta; zetai++)
    {
        canonicalMoments_[zetai]
            = zeta[zetai]/(1.0 - canonicalMoments_[zetai - 1]);

        if 
        (
            canonicalMoments_[zetai] < smallZeta_ 
         || canonicalMoments_[zetai] > 1.0
        )
        {
            nRealizableMoments_ = zetai + 1;

            return;
        }
        else if
        (
            mag(canonicalMoments_[zetai]) <= smallZeta_
         || mag(canonicalMoments_[zetai] - 1.0) <= smallZeta_
        )
        {
            nRealizableMoments_ = zetai + 2;
            onMomentSpaceBoundary_ = true;

            return;
        }
    }

    onMomentSpaceBoundary_ = false;
    nRealizableMoments_ = nZeta + 1;
}

void Foam::univariateMomentSet::checkRealizability
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (realizabilityChecked_)
    {
        return;
    }

    // If the zero-order moment is negative, exit immediately.
    if ((*this)[0] < 0.0)
    {
        if (fatalErrorOnFailedRealizabilityTest)
        {
            // If the user requested to throw an error when the realizability
            // test fails, we do so.
            FatalErrorInFunction
                << "The zero-order moment is negative." << nl
                << "    Moment set: " << (*this)
                << abort(FatalError);
        }
        else
        {
            // If the zero-order moment is negative, the moment set is not 
            // realizable. If the user has requested to not throw an error,
            // we mark the moment set as not realizable, set the number of 
            // realizable moments to zero and return. This is necessary when
            // using some adaptive methods which explicitly test for
            // realizability to make decisions.
            realizabilityChecked_ = true;
            negativeZeta_ = 0;
            nRealizableMoments_ = 0;
            fullyRealizable_ = false;
            subsetRealizable_ = false;
            onMomentSpaceBoundary_ = false;

            return;
        }
    }

    // Set flags and return if the zero-order moment is too small but an
    // error should not be thrown. Do nothing otherwise.
    if ((*this)[0] < smallM0_ && !fatalErrorOnFailedRealizabilityTest)
    {
        realizabilityChecked_ = true;
        negativeZeta_ = 0;
        nRealizableMoments_ = 0;
        fullyRealizable_ = false;
        subsetRealizable_ = false;
        onMomentSpaceBoundary_ = false;

        return;
    }

    // Check for the degenerate case where only m0 is defined and throw an error
    // if this is the case.
    if (nMoments_ <= 1)
    {
        FatalErrorInFunction
            << "The moment has size less or equal to 1." << nl
            << "    Moment set: " << (*this)
            << abort(FatalError);
    }

    // Reset vector of zeta values to check realizability
    zeta_ = 0.0;

    // Check for the case with only two moments
    if (nMoments_ == 2)
    {
        // In the case of support over R, m0 must be positive and m1 needs to
        // be a real number. The first condition is satisfied, so only flags
        // need to be set before returning.
        if (support_ == "R")
        {
            realizabilityChecked_ = true;
            negativeZeta_ = 0;
            nRealizableMoments_ = 2;
            fullyRealizable_ = true;
            subsetRealizable_ = true;
            onMomentSpaceBoundary_ = false;

            return;
        }

        // Managing other supports (R+ and [0, 1])

        // Calculate zeta_1 (we do not store zeta_0 because it is always 1)
        zeta_[0] = (*this)[1]/(*this)[0];

        if (zeta_[0] <= smallZeta_)
        {
            if (isDegenerate() || zeta_[0] == 0.0)
            {
                realizabilityChecked_ = true;
                negativeZeta_ = 0;
                nRealizableMoments_ = 2;
                fullyRealizable_ = true;
                subsetRealizable_ = true;
                onMomentSpaceBoundary_ = true;

                return;
            }

            if (fatalErrorOnFailedRealizabilityTest)
            {
                FatalErrorInFunction
                << "Moment set with dimension 2 and only one valid moment."
                << nl << "    Moment set: " << (*this)
                << abort(FatalError);
            }
            else
            {
                realizabilityChecked_ = true;
                negativeZeta_ = 1;
                nRealizableMoments_ = 1;
                fullyRealizable_ = false;
                subsetRealizable_ = false;
                onMomentSpaceBoundary_ = false;

                return;
            }
        }

        if (support_ == "RPlus") // Support on R+ - Check if zetas are positive.
        {
            realizabilityChecked_ = true;
            negativeZeta_ = 0;
            nRealizableMoments_ = 2;
            fullyRealizable_ = true;
            subsetRealizable_ = true;
            onMomentSpaceBoundary_ = false;

            return;
        }
        else // Support on [0, 1] - Check if canonical moments belong to [0,1].
        {
            if (zeta_[0] <= 1.0)
            {
                realizabilityChecked_ = true;
                nRealizableMoments_ = 2;
                fullyRealizable_ = true;
                subsetRealizable_ = true;

                if (zeta_[0] < 1.0)
                {
                    onMomentSpaceBoundary_ = false;
                }
                else
                {
                    onMomentSpaceBoundary_ = true;
                }

                return;
            }
            else
            {
                if (isDegenerate())
                {
                    realizabilityChecked_ = true;
                    negativeZeta_ = 0;
                    nRealizableMoments_ = 2;
                    fullyRealizable_ = true;
                    subsetRealizable_ = true;
                    onMomentSpaceBoundary_ = true;

                    return;
                }

                if (fatalErrorOnFailedRealizabilityTest)
                {
                    FatalErrorInFunction
                    << "Moment set with dimension 2 and only one valid moment."
                    << nl << "    Moment set: " << (*this)
                    << abort(FatalError);
                }
                else
                {
                    realizabilityChecked_ = true;
                    negativeZeta_ = 1;
                    nRealizableMoments_ = 1;
                    fullyRealizable_ = false;
                    subsetRealizable_ = false;
                    onMomentSpaceBoundary_ = false;

                    return;
                }
            }
        }
    }

    // Check for the case with more than two moments

    // Store the number of zeta elements.
    // This is the number of moments minus one, but it was already calculated
    // in the constructor, when zeta_ was initialized. It is copied to ensure
    // consistency.
    label nN = zeta_.size();

    // Calculate the integer part of nN/2.
    label nD = label(nN/2);

    // Calculate the remainder of the division nN/2.
    label nR = nN - 2*nD;

    // Matrix used to build the recurrence relation
    scalarRectangularMatrix zRecurrence(nD + 1, nMoments_, Zero);

    for (label columnI = 0; columnI < nMoments_; columnI++)
    {
        zRecurrence[0][columnI] = (*this)[columnI]/(*this)[0];
    }

    alpha_[0] = (*this)[1]/(*this)[0];
    beta_[0] = 1.0;

    for (label columnI = 1; columnI < nMoments_ - 1; columnI++)
    {
        zRecurrence[1][columnI] = zRecurrence[0][columnI + 1]
              - alpha_[0]*zRecurrence[0][columnI];
    }

    zeta_[0] = alpha_[0];

    if (!(support_ == "R") && zeta_[0] <= smallZeta_)
    {
        if (isDegenerate() || zeta_[0] == 0.0)
        {
            realizabilityChecked_ = true;
            negativeZeta_ = 0;
            nRealizableMoments_ = 2;
            fullyRealizable_ = false;
            subsetRealizable_ = true;
            onMomentSpaceBoundary_ = true;

            return;
        }

        if (fatalErrorOnFailedRealizabilityTest)
        {
            FatalErrorInFunction
                << "Moment set with only one valid moment."
                << nl << "    Moment set: " << (*this)
                << "zeta vector = " << zeta_ << endl
                << "smallZeta = " << smallZeta_
                << abort(FatalError);
        }
        else
        {
            realizabilityChecked_ = true;
            negativeZeta_ = 1;
            nRealizableMoments_ = 1;
            fullyRealizable_ = false;
            subsetRealizable_ = false;
            onMomentSpaceBoundary_ = false;

            return;
        }
    }

    for (label zetai = 1; zetai < nD; zetai++)
    {
        beta_[zetai] = zRecurrence[zetai][zetai]
                /zRecurrence[zetai - 1][zetai - 1];

        if (support_ == "R")
        {
            if (beta_[zetai] <= smallZeta_)
            {
                realizabilityChecked_ = true;
                nRealizableMoments_ = 2*zetai;
                fullyRealizable_ = false;
                subsetRealizable_ = true;

                return;
            }
        }
        else
        {
            zeta_[2*zetai - 1] = beta_[zetai]/zeta_[2*zetai - 2];

            if (zeta_[2*zetai - 1] <= smallZeta_)
            {
                if (support_ == "RPlus")
                {
                    if (zeta_[2*zetai - 1] < smallZeta_)
                    {
                        negativeZeta_ = 2*zetai;
                        nRealizableMoments_ = negativeZeta_;
                        onMomentSpaceBoundary_ = false;
                    }
                    else
                    {
                        negativeZeta_ = 2*zetai + 1;
                        nRealizableMoments_ = negativeZeta_;
                        onMomentSpaceBoundary_ = true;
                    }
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta_, 2*zetai);
                }

                realizabilityChecked_ = true;
                fullyRealizable_ = false;
                subsetRealizable_ = true;

                return;
            }
        }

        alpha_[zetai] =
            zRecurrence[zetai][zetai + 1]
           /max(zRecurrence[zetai][zetai], SMALL)
          - zRecurrence[zetai - 1][zetai]
           /zRecurrence[zetai - 1][zetai - 1];

        if (!(support_ == "R"))
        {
            zeta_[2*zetai] = alpha_[zetai] - zeta_[2*zetai - 1];

            if (zeta_[2*zetai] <= smallZeta_)
            {
                if (support_ == "RPlus")
                {
                    if (zeta_[2*zetai] < smallZeta_)
                    {
                        negativeZeta_ = 2*zetai + 1;
                        nRealizableMoments_ = negativeZeta_;
                        onMomentSpaceBoundary_ = false;
                    }
                    else
                    {
                        negativeZeta_ = 2*zetai + 2;
                        nRealizableMoments_ = negativeZeta_;
                        onMomentSpaceBoundary_ = true;
                    }
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta_, 2*zetai + 1);
                }

                realizabilityChecked_ = true;
                fullyRealizable_ = false;
                subsetRealizable_ = true;

                return;
            }
        }

        for (label columnI = zetai + 1; columnI <= nN - zetai - 1; columnI++)
        {
            zRecurrence[zetai + 1][columnI] = zRecurrence[zetai][columnI + 1]
                    - alpha_[zetai]*zRecurrence[zetai][columnI]
                    - beta_[zetai]*zRecurrence[zetai - 1][columnI];
        }
    }
   
    beta_[nD] = zRecurrence[nD][nD]/max(zRecurrence[nD - 1][nD - 1], SMALL);

    if (support_ == "R")
    {
        alpha_[nD] = zRecurrence[nD][nD + 1]/max(zRecurrence[nD][nD], SMALL)
                    - zRecurrence[nD - 1][nD]/max(zRecurrence[nD - 1][nD - 1], 
                    SMALL);

        if (beta_[nD] <= smallZeta_)
        {
            realizabilityChecked_ = true;
            nRealizableMoments_ = 2*nD;
            fullyRealizable_ = false;
            subsetRealizable_ = true;

            return;
        }
        else
        {
            realizabilityChecked_ = true;
            nRealizableMoments_ = nMoments_;
            fullyRealizable_ = true;
            subsetRealizable_ = true;

            return;
        }
    }
    else
    {
        zeta_[2*nD - 1] = beta_[nD]/zeta_[2*nD - 2];

        if (zeta_[2*nD - 1] <= smallZeta_)
        {
            if (support_ == "RPlus")
            {
                if (zeta_[2*nD - 1] < smallZeta_)
                {
                    negativeZeta_ = 2*nD;
                    nRealizableMoments_ = negativeZeta_;
                    onMomentSpaceBoundary_ = false;
                }
                else
                {
                    negativeZeta_ = 2*nD + 1;
                    nRealizableMoments_ = negativeZeta_;
                    onMomentSpaceBoundary_ = true;
                }
            }
            else  // Support on [0,1]
            {
                checkCanonicalMoments(zeta_, 2*nD);
            }

            realizabilityChecked_ = true;
            fullyRealizable_ = false;
            subsetRealizable_ = true;

            return;
        }

        if (nR == 1)
        {
            alpha_[nD] = zRecurrence[nD][nD + 1]/zRecurrence[nD][nD]
                    - zRecurrence[nD - 1][nD]/zRecurrence[nD - 1][nD - 1];

            zeta_[2*nD] = alpha_[nD] - zeta_[2*nD - 1];

            if (zeta_[2*nD] <= smallZeta_)
            {
                if (support_ == "RPlus")
                {
                    if (zeta_[2*nD] < smallZeta_)
                    {
                        negativeZeta_ = 2*nD + 1;
                        nRealizableMoments_ = negativeZeta_;
                        fullyRealizable_ = false;
                        onMomentSpaceBoundary_ = false;
                    }
                    else
                    {
                        negativeZeta_ = nN;
                        nRealizableMoments_ = nMoments_;
                        fullyRealizable_ = true;
                        onMomentSpaceBoundary_ = true;
                    }
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta_, 2*nD + 1);

                    if (onMomentSpaceBoundary_)
                    {
                        negativeZeta_ = nN;
                        fullyRealizable_ = true;
                    }
                    else
                    {
                        negativeZeta_ = 2*nD + 1;
                        fullyRealizable_ = false;
                    }
                }

                realizabilityChecked_ = true;
                subsetRealizable_ = true;

                return;
            }
            else // zeta_[2*nD] > 0.0
            {
                if (support_ == "RPlus")
                {
                    negativeZeta_ = nN;
                    nRealizableMoments_ = nMoments_;
                    fullyRealizable_ = true;
                    onMomentSpaceBoundary_ = false;
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta_, 2*nD + 1);

                    if (onMomentSpaceBoundary_)
                    {
                        negativeZeta_ = nN;
                        fullyRealizable_ = true;
                    }
                    else
                    {
                        negativeZeta_ = 2*nD + 1;
                        fullyRealizable_ = false;
                    }
                }

                realizabilityChecked_ = true;
                subsetRealizable_ = true;

                return;
            }
        }
        else
        {
            // If support is [0, + inf[ and this level is reached, the full set
            // of moments is realizable
            if (support_ == "RPlus")
            {
                negativeZeta_ = nN;
                fullyRealizable_ = true;
                subsetRealizable_ = true;
                nRealizableMoments_ = nMoments_;
                onMomentSpaceBoundary_ = false;
            }
            else
            {
                checkCanonicalMoments(zeta_, nN);

                if (nRealizableMoments_ == nMoments_)
                {
                    negativeZeta_ = nN;
                    fullyRealizable_ = true;
                    subsetRealizable_ = true;
                }
                else
                {
                    negativeZeta_ = nN;
                    fullyRealizable_ = false;
                    subsetRealizable_ = true;
                }
            }

            realizabilityChecked_ = true;

            return;
        }
    }
}


Foam::labelListList Foam::univariateMomentSet::makeUnivariateMomentOrders
(
    const label nMoments
)
{
    labelListList mOrd(nMoments);

    for (label mI = 0; mI < nMoments; mI++)
    {
        mOrd[mI] = labelList(1, mI);
    }

    return mOrd;
}


void Foam::univariateMomentSet::update
(
    const scalarList& weights,
    const scalarList& abscissae
)
{
    // Recomputing all the moments (even if they originally were not realizable)
    // from quadrature (projection step).
    for (label momenti = 0; momenti < nMoments_; momenti++)
    {
        (*this)[momenti] = Zero;

        for (label nodei = 0; nodei < weights.size(); nodei++)
        {
            (*this)[momenti] += weights[nodei]*pow(abscissae[nodei], momenti);
        }
    }

    realizabilityChecked_ = false;
}


void Foam::univariateMomentSet::setSize(const label newSize)
{
    label oldSize = (*this).size();
    Foam::momentSet::setSize(newSize);
    realizabilityChecked_ = false;

    if (oldSize > newSize)
    {
        makeUnivariateMomentOrders(newSize);
    }
}


void Foam::univariateMomentSet::resize(const label newSize)
{
    setSize(newSize);
}

// ************************************************************************* //
