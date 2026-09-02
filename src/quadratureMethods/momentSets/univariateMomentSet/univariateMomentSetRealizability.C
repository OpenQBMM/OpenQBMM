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

Description
    Realizability check of a univariate moment set, and construction of the
    recurrence relationship it determines.

\*---------------------------------------------------------------------------*/

#include "univariateMomentSet.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::univariateMomentSet::checkCanonicalMoments
(
    const scalarList& zeta,
    const label nZeta
)
{
    canonicalMoments_[0] = zeta[0];

    if (mag(canonicalMoments_[0] - 1.0) <= smallZeta_)
    {
        status_.nRealizableMoments = 2;
        status_.onMomentSpaceBoundary = true;

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
            status_.nRealizableMoments = zetai + 1;

            return;
        }
        else if
        (
            mag(canonicalMoments_[zetai]) <= smallZeta_
         || mag(canonicalMoments_[zetai] - 1.0) <= smallZeta_
        )
        {
            status_.nRealizableMoments = zetai + 2;
            status_.onMomentSpaceBoundary = true;

            return;
        }
    }

    status_.onMomentSpaceBoundary = false;
    status_.nRealizableMoments = nZeta + 1;
}


bool Foam::univariateMomentSet::checkTrivialCases
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    const supportType& mSupport = supports_[0];
    const scalar m0 = moment(0);

    // If the zero-order moment is negative, the moment set is not realizable.
    if (m0 < 0.0)
    {
        if (fatalErrorOnFailedRealizabilityTest)
        {
            FatalErrorInFunction
                << "The zero-order moment is negative." << nl
                << "    Moment set: " << *this << nl
                << exit(FatalError);
        }

        // If the user has requested to not throw an error, the moment set is
        // marked as not realizable and the number of realizable moments is
        // set to zero. This is necessary when using adaptive methods which
        // explicitly test for realizability to make decisions.
        status_.checked = true;
        status_.fullyRealizable = false;
        status_.subsetRealizable = false;

        return true;
    }

    // A negligible zero-order moment carries no information. An error is
    // thrown only if the caller asked for it, so that adaptive methods can
    // skip the moment set instead.
    if (m0 < smallM0_ && !fatalErrorOnFailedRealizabilityTest)
    {
        status_.checked = true;
        status_.fullyRealizable = false;
        status_.subsetRealizable = false;

        return true;
    }

    // Check for the degenerate case where only m0 is defined
    if (nMoments() <= 1)
    {
        FatalErrorInFunction
            << "The moment set has size less or equal to 1." << nl
            << "    Moment set: " << *this << nl
            << exit(FatalError);
    }

    if (nMoments() > 2)
    {
        return false;
    }

    // Moment set with two moments. The single-node quadrature it defines is
    // fully determined, so the first coefficients of the recurrence
    // relationship are set here as well.
    const scalar zeta0 = moment(1)/m0;

    alpha_[0] = zeta0;
    beta_[0] = 1.0;
    zeta_[0] = zeta0;

    if (mSupport == supportType::ZeroOne)
    {
        canonicalMoments_[0] = zeta0;
    }

    status_.checked = true;
    status_.subsetRealizable = true;

    // Over R, m0 must be positive and m1 needs to be a real number. Both
    // conditions are satisfied at this point.
    if (mSupport == supportType::R)
    {
        status_.nRealizableMoments = 2;
        status_.fullyRealizable = true;

        return true;
    }

    if (zeta0 <= smallZeta_)
    {
        if (isDegenerate() || zeta0 == 0.0)
        {
            status_.nRealizableMoments = 2;
            status_.fullyRealizable = true;
            status_.onMomentSpaceBoundary = true;

            return true;
        }

        if (fatalErrorOnFailedRealizabilityTest)
        {
            FatalErrorInFunction
                << "Moment set with dimension 2 and only one valid moment."
                << nl
                << "    Moment set: " << *this << nl
                << exit(FatalError);
        }

        status_.nRealizableMoments = 1;
        status_.fullyRealizable = false;
        status_.subsetRealizable = false;

        return true;
    }

    // Support over R+: the moment set is realizable because zeta_0 > 0
    if (mSupport == supportType::RPlus)
    {
        status_.nRealizableMoments = 2;
        status_.fullyRealizable = true;

        return true;
    }

    // Support over [0, 1]: the canonical moment must belong to [0, 1]
    if (zeta0 <= 1.0)
    {
        status_.nRealizableMoments = 2;
        status_.fullyRealizable = true;
        status_.onMomentSpaceBoundary = (zeta0 == 1.0);

        return true;
    }

    if (isDegenerate())
    {
        status_.nRealizableMoments = 2;
        status_.fullyRealizable = true;
        status_.onMomentSpaceBoundary = true;

        return true;
    }

    if (fatalErrorOnFailedRealizabilityTest)
    {
        FatalErrorInFunction
            << "Moment set with dimension 2 and only one valid moment." << nl
            << "    Moment set: " << *this << nl
            << exit(FatalError);
    }

    status_.nRealizableMoments = 1;
    status_.fullyRealizable = false;
    status_.subsetRealizable = false;

    return true;
}


void Foam::univariateMomentSet::initialiseRecurrence()
{
    const label nM = nMoments();
    const label nD = nZeta()/2;
    const scalar m0 = moment(0);

    // Resize the Wheeler table only if necessary, so that a moment set reused
    // over the cells of a mesh does not reallocate it
    if (zRecurrence_.m() != nD + 1 || zRecurrence_.n() != nM)
    {
        zRecurrence_.setSize(nD + 1, nM);
    }

    zRecurrence_ = 0.0;

    for (label columnI = 0; columnI < nM; columnI++)
    {
        zRecurrence_[0][columnI] = moment(columnI)/m0;
    }

    alpha_[0] = moment(1)/m0;
    beta_[0] = 1.0;
    zeta_[0] = alpha_[0];

    for (label columnI = 1; columnI < nM - 1; columnI++)
    {
        zRecurrence_[1][columnI] =
            zRecurrence_[0][columnI + 1] - alpha_[0]*zRecurrence_[0][columnI];
    }
}


bool Foam::univariateMomentSet::checkFirstZeta
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    // The zeta_k are not computed for measures with support over R
    if (supports_[0] == supportType::R || zeta_[0] > smallZeta_)
    {
        return false;
    }

    status_.checked = true;

    if (isDegenerate() || zeta_[0] == 0.0)
    {
        status_.nRealizableMoments = 2;
        status_.fullyRealizable = false;
        status_.subsetRealizable = true;
        status_.onMomentSpaceBoundary = true;

        return true;
    }

    if (fatalErrorOnFailedRealizabilityTest)
    {
        FatalErrorInFunction
            << "Moment set with only one valid moment." << nl
            << "    Moment set: " << *this << nl
            << "    Zeta vector: " << zeta_ << nl
            << "    smallZeta: " << smallZeta_ << nl
            << exit(FatalError);
    }

    status_.nRealizableMoments = 1;
    status_.fullyRealizable = false;
    status_.subsetRealizable = false;

    return true;
}


Foam::label Foam::univariateMomentSet::buildRecurrence()
{
    const bool overR = (supports_[0] == supportType::R);
    const label nZ = nZeta();
    const label nD = nZ/2;

    for (label j = 1; j <= nD; j++)
    {
        // The denominators of the Wheeler recursion are strictly positive
        // wherever the walk reaches them: the table is seeded with
        // zRecurrence_[0][0] = 1, and the positivity the previous step
        // required of beta_{j-1}, directly for support over R and through
        // zeta_{2j-3} otherwise, is what makes zRecurrence_[j-1][j-1]
        // positive. They are therefore divided by as they are. Flooring them
        // would leave the quotient of a moment set that is close to
        // degenerate, but still realizable, silently wrong.
        beta_[j] = zRecurrence_[j][j]/zRecurrence_[j - 1][j - 1];

        // Odd position of the zeta chain, determined by beta_j
        const label oddZetai = 2*j - 1;

        if (overR)
        {
            // Over R only the positivity of the beta_k is required
            if (beta_[j] <= smallZeta_)
            {
                return oddZetai;
            }
        }
        else
        {
            zeta_[oddZetai] = beta_[j]/zeta_[oddZetai - 1];

            if (zeta_[oddZetai] <= smallZeta_)
            {
                return oddZetai;
            }
        }

        // Even position of the zeta chain, determined by alpha_j. It exists
        // only if the moments determine alpha_j: with an odd number of
        // moments, alpha_nD would require the moment of order nMoments, which
        // is not available.
        const label evenZetai = 2*j;

        if (evenZetai <= nZ - 1)
        {
            alpha_[j] =
                zRecurrence_[j][j + 1]/zRecurrence_[j][j]
              - zRecurrence_[j - 1][j]/zRecurrence_[j - 1][j - 1];

            if (!overR)
            {
                zeta_[evenZetai] = alpha_[j] - zeta_[oddZetai];

                if (zeta_[evenZetai] <= smallZeta_)
                {
                    return evenZetai;
                }
            }
        }

        // Advance the Wheeler table to the next row
        if (j < nD)
        {
            for (label columnI = j + 1; columnI <= nZ - j - 1; columnI++)
            {
                zRecurrence_[j + 1][columnI] =
                    zRecurrence_[j][columnI + 1]
                  - alpha_[j]*zRecurrence_[j][columnI]
                  - beta_[j]*zRecurrence_[j - 1][columnI];
            }
        }
    }

    return nZ;
}


void Foam::univariateMomentSet::classify
(
    const label zetai,
    const bool complete
)
{
    const supportType& mSupport = supports_[0];
    const label nZ = nZeta();

    status_.checked = true;
    status_.subsetRealizable = true;

    if (mSupport == supportType::ZeroOne)
    {
        // The canonical moments determine how many moments are realizable and
        // whether the moment set lies on the boundary of the moment space
        checkCanonicalMoments(zeta_, complete ? nZ : zetai + 1);
    }
    else if (complete)
    {
        status_.nRealizableMoments = nMoments();
    }
    else if (mSupport == supportType::R)
    {
        status_.nRealizableMoments = zetai + 1;
    }
    else
    {
        // A zeta_k that is null within the tolerance leaves the moment set on
        // the boundary of the moment space, with one more realizable moment
        // than a strictly negative one
        status_.onMomentSpaceBoundary = !(zeta_[zetai] < smallZeta_);

        status_.nRealizableMoments =
            zetai + (status_.onMomentSpaceBoundary ? 2 : 1);
    }

    status_.fullyRealizable = (status_.nRealizableMoments == nMoments());
}


void Foam::univariateMomentSet::checkRealizability
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (status_.checked)
    {
        return;
    }

    // Reset the quantities derived from the moment vector, so that the
    // outcome does not depend on the previous contents of a moment set that
    // is reused
    resetStatus();

    if (checkTrivialCases(fatalErrorOnFailedRealizabilityTest))
    {
        return;
    }

    initialiseRecurrence();

    if (checkFirstZeta(fatalErrorOnFailedRealizabilityTest))
    {
        return;
    }

    const label zetai = buildRecurrence();

    classify(zetai, zetai == nZeta());
}


// ************************************************************************* //
