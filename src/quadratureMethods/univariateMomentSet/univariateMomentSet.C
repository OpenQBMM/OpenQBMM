/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 Alberto Passalacqua
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

#include "univariateMomentSet.H"

#include "eigenSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMomentSet::univariateMomentSet
(
    const label nMoments,
    const scalar initValue,
    const word support
)
:
    scalarDiagonalMatrix(nMoments, initValue),
    nMoments_(nMoments),
    alpha_(label((nMoments_ - 2)/2) + 1, scalar(0)),
    beta_(label((nMoments_ - 1)/2) + 1, scalar(0)),
    support_(support),
    inverted_(false),
    fullyRealizable_(true),
    subsetRealizable_(true),
    realizabilityChecked_(false),
    quadratureSetUp_(false),
    nInvertibleMoments_(nMoments_)
{
    if (support_ != "R" && support_ != "RPlus" && support_ != "01")
    {
        FatalErrorIn
        (
            "Foam::univariateMomentSet::univariateMomentSet\n"
            "(\n"
            "    const label nMoments,\n"
            "    const scalar initValue\n"
            "    const word support\n"
            ")"
        )   << "The specified support is invalid.\n"
            << "Valid supports are: R, RPlus and 01."
            << abort(FatalError);
    }
}

Foam::univariateMomentSet::univariateMomentSet
(
    const scalarDiagonalMatrix& m,
    const word support
)
:
    scalarDiagonalMatrix(m),
    nMoments_(m.size()),
    alpha_(label((nMoments_ - 2)/2) + 1, scalar(0)),
    beta_(label((nMoments_ - 1)/2) + 1, scalar(0)),
    support_(support),
    inverted_(false),
    fullyRealizable_(true),
    subsetRealizable_(true),
    realizabilityChecked_(false),
    quadratureSetUp_(false),
    nInvertibleMoments_(nMoments_)
{
    if (support_ != "R" && support_ != "RPlus" && support_ != "01")
    {
        FatalErrorIn
        (
            "Foam::univariateMomentSet::univariateMomentSet\n"
            "(\n"
            "    const scalarDiagonalMatrix& m\n"
            "    const word support\n"
            ")"
        )   << "The specified support is invalid.\n"
            << "Valid supports are: R, RPlus and 01."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateMomentSet::~univariateMomentSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::univariateMomentSet::invert()
{
    if (inverted_)
    {
        return;
    }

    if (!realizabilityChecked_)
    {
        checkRealizability();
        setupQuadrature(true);
    }

    if (!quadratureSetUp_)
    {
        setupQuadrature(true);
    }

    if (nInvertibleMoments_ < 2)
    {
        FatalErrorIn
        (
            "Foam::univariateMomentSet::invert\n"
            "(\n"
            "    const scalarDiagonalMatrix& weights,\n"
            "    const scalarDiagonalMatrix& abscissae\n"
            ")"
        )   << "Insufficient number (" << nInvertibleMoments_
            << ") of moments to define quadrature."
            << abort(FatalError);
    }

    if (nInvertibleMoments_ == 2)
    {
        weights_[0] = 1.0;
        abscissae_[0] = (*this)[1]/(*this)[0];

        inverted_ = true;

        return;
    }

    scalarSquareMatrix z(nNodes_, nNodes_, scalar(0));

    for (label i = 0; i < nNodes_ - 1; i++)
    {
        z[i][i] = alpha_[i];
        z[i][i+1] = Foam::sqrt(beta_[i+1]);
        z[i+1][i] = z[i][i+1];
    }

    z[nNodes_ - 1][nNodes_ - 1] = alpha_[nNodes_ - 1];

    // Computing weights and abscissae
    eigenSolver zEig(z, true);

    // Computing weights and abscissae
    for (label i = 0; i < nNodes_; i++)
    {
        weights_[i] = (*this)[0]*sqr(zEig.eigenvectors()[0][i]);
        abscissae_[i] = zEig.eigenvaluesRe()[i];
    }

    inverted_ = true;
}

void Foam::univariateMomentSet::checkCanonicalMoments
(
    const scalarDiagonalMatrix& zeta,
    const label nZeta,
    label& nRealizableMoments
)
{
    Info << endl << "Checking canonical moments" << endl;
    Info << "nZeta = " << nZeta << endl;

    scalarDiagonalMatrix canonicalMoments(nZeta, 0.0);

    canonicalMoments[0] = zeta[0];

    for (label zetaI = 1; zetaI < nZeta; zetaI++)
    {
        Info << "Loop " << zetaI << endl;

        canonicalMoments[zetaI] = zeta[zetaI]/(1.0 - canonicalMoments[zetaI]);

        if (canonicalMoments[zetaI] < 0 || canonicalMoments[zetaI] > 1)
        {
            nRealizableMoments = zetaI + 1;

            return;
        }
    }

    nRealizableMoments = nZeta + 1;

    Info << "Realizable moments: " << nRealizableMoments << endl;
}

void Foam::univariateMomentSet::checkRealizability()
{
    if (realizabilityChecked_)
    {
        return;
    }

    // If the zero-order moment is negative, exit immediately.
    if ((*this)[0] < 0)
    {
        FatalErrorIn
        (
            "Foam::univariateMomentSet::checkRealizability()\n"
        )   << "The zero-order moment is negative."
            << abort(FatalError);
    }

    // Check for the degenerate case where only m0 is defined
    if (nMoments_ <= 1)
    {
        FatalErrorIn
        (
            "Foam::univariateMomentSet::checkRealizability()\n"
        )   << "The moment set is degenerate."
            << abort(FatalError);
    }

    label nN = nMoments_ - 1;
    label nD = label(nN/2);
    label nR = nN - 2*nD;

    // Vector of zeta values used to check moment realizability
    scalarDiagonalMatrix zeta(nN);

    // Check for the case with only two moments, if support is R+ or [0,1]
    // In the case of support over R, only check if beta < 0
    if (nMoments_ == 2)
    {
        if (support_ == "R")
        {
            negativeZeta_ = 0;
            nRealizableMoments_ = 2;
            nInvertibleMoments_ = 2;
            fullyRealizable_ = true;
            subsetRealizable_ = true;

            return;
        }

        zeta[0] = (*this)[1]/(*this)[0];

        if (zeta[0] <= 0)
        {
            negativeZeta_ = 1;
            nRealizableMoments_ = 1;
            nInvertibleMoments_ = 0;
            fullyRealizable_ = false;
            subsetRealizable_ = true;

            FatalErrorIn
            (
                "Foam::univariateMomentSet::checkRealizability()\n"
            )   << "Moment set with dimension 2 and only one realizable moment."
                << abort(FatalError);
        }

        if (support_ == "RPlus")
        {
            negativeZeta_ = 0;
            nRealizableMoments_ = 2;
            nInvertibleMoments_ = 2;
            fullyRealizable_ = true;
            subsetRealizable_ = true;
            realizabilityChecked_ = true;

            return;
        }
        else // Support on [0, 1] - Check if canonical moments belong to [0,1]
        {
            if (zeta[0] <= 1.0)
            {
                nRealizableMoments_ = 2;
                nInvertibleMoments_ = 2;
                fullyRealizable_ = true;
                subsetRealizable_ = true;
                realizabilityChecked_ = true;

                return;
            }
            else
            {
                negativeZeta_ = 1;
                nRealizableMoments_ = 1;
                nInvertibleMoments_ = 0;
                fullyRealizable_ = false;
                subsetRealizable_ = true;

                FatalErrorIn
                (
                    "Foam::univariateMomentSet::checkRealizability()\n"
                )   << "Moment set with dimension 2 and only one "
                    << "realizable moment."
                    << abort(FatalError);
            }
        }
    }

    // Matrix used to build the recurrence relation
    scalarRectangularMatrix zRecurrence(nD + 1, nMoments_, 0.0);

    for (label columnI = 0; columnI < nMoments_; columnI++)
    {
        zRecurrence[0][columnI] = (*this)[columnI];
    }

    alpha_[0] = (*this)[1]/(*this)[0];
    beta_[0] = (*this)[0];

    for (label columnI = 1; columnI < nMoments_ - 1; columnI++)
    {
        zRecurrence[1][columnI] = zRecurrence[0][columnI + 1]
              - alpha_[0]*zRecurrence[0][columnI];
    }

    zeta[0] = alpha_[0];

    if (!(support_ == "R") && zeta[0] <= 0)
    {
        negativeZeta_ = 1;
        nRealizableMoments_ = 1;
        nInvertibleMoments_ = 0;
        fullyRealizable_ = false;
        subsetRealizable_ = true;

        FatalErrorIn
        (
            "Foam::univariateMomentSet::checkRealizability()\n"
        )   << "Moment set with only one realizable moment."
            << abort(FatalError);

        // No need to check canonical moments
    }

    for (label zetaI = 1; zetaI <= nD - 1; zetaI++)
    {
        beta_[zetaI] = zRecurrence[zetaI][zetaI]
                /zRecurrence[zetaI - 1][zetaI - 1];

        if (support_ == "R" && beta_[zetaI] < 0.0)
        {
            nRealizableMoments_ = 2*zetaI;
            calcNInvertibleMoments();
            fullyRealizable_ = false;
            subsetRealizable_ = true;
            realizabilityChecked_ = true;

            return;
        }
        else
        {
            zeta[2*zetaI - 1] = beta_[zetaI]/zeta[2*zetaI - 2];

            if (zeta[2*zetaI - 1] <= 0)
            {
                if (support_ == "RPlus")
                {
                    negativeZeta_ = 2*zetaI;
                    nRealizableMoments_ = negativeZeta_;
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta, 2*zetaI - 1, nRealizableMoments_);
                }

                calcNInvertibleMoments();
                fullyRealizable_ = false;
                subsetRealizable_ = true;
                realizabilityChecked_ = true;

                return;
            }
        }

        alpha_[zetaI] = zRecurrence[zetaI][zetaI + 1]/zRecurrence[zetaI][zetaI]
                - zRecurrence[zetaI - 1][zetaI]
                /zRecurrence[zetaI - 1][zetaI - 1];

        if (!(support_ == "R"))
        {
            zeta[2*zetaI] = alpha_[zetaI] - zeta[2*zetaI - 1];

            if (zeta[2*zetaI] <= 0)
            {
                if (support_ == "RPlus")
                {
                    negativeZeta_ = 2*zetaI + 1;
                    nRealizableMoments_ = negativeZeta_;
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta, 2*zetaI, nRealizableMoments_);
                }

                calcNInvertibleMoments();
                fullyRealizable_ = false;
                subsetRealizable_ = true;
                realizabilityChecked_ = true;

                return;
            }
        }

        for (label columnI = zetaI + 1; columnI <= nN - zetaI - 1; columnI++)
        {
            zRecurrence[zetaI + 1][columnI] = zRecurrence[zetaI][columnI + 1]
                    - alpha_[zetaI]*zRecurrence[zetaI][columnI]
                    - beta_[zetaI]*zRecurrence[zetaI - 1][columnI];
        }
    }

    beta_[nD] = zRecurrence[nD][nD]/zRecurrence[nD - 1][nD - 1];

    if (support_ == "R")
    {
        if (beta_[nD] < 0.0)
        {
            nRealizableMoments_ = 2*nD;
            calcNInvertibleMoments();
            fullyRealizable_ = false;
            subsetRealizable_ = true;
            realizabilityChecked_ = true;

            return;
        }
        else
        {
            nRealizableMoments_ = nMoments_;
            calcNInvertibleMoments();
            fullyRealizable_ = true;
            subsetRealizable_ = true;
            realizabilityChecked_ = true;

            return;
        }
    }
    else
    {
        zeta[2*nD - 1] = beta_[nD]/zeta[2*nD - 2];

        if (zeta[2*nD - 1] <= 0)
        {
            Info << "Checking CM OUT 2dn - 1" << endl;
            if (support_ == "RPlus")
            {
                negativeZeta_ = 2*nD;
                nRealizableMoments_ = negativeZeta_;
            }
            else  // Support on [0,1]
            {
                checkCanonicalMoments(zeta, 2*nD - 1, nRealizableMoments_);
            }

            calcNInvertibleMoments();
            fullyRealizable_ = false;
            subsetRealizable_ = true;
            realizabilityChecked_ = true;

            return;
        }

        if (nR == 1)
        {
            alpha_[nD] = zRecurrence[nD][nD + 1]/zRecurrence[nD][nD]
                    - zRecurrence[nD - 1][nD]/zRecurrence[nD - 1][nD - 1];

            zeta[2*nD] = alpha_[nD] - zeta[2*nD - 1];

            if (zeta[2*nD] <= 0)
            {
                if (support_ == "RPlus")
                {
                    negativeZeta_ = 2*nD + 1;
                    nRealizableMoments_ = negativeZeta_;
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta, 2*nD, nRealizableMoments_);
                }

                calcNInvertibleMoments();
                fullyRealizable_ = false;
                subsetRealizable_ = true;
                realizabilityChecked_ = true;

                return;
            }
        }
        else
        {
            // If support is [0, + inf[ and this level is reached, the full set of
            // moments is realizable
            if (support_ == "RPlus")
            {
                fullyRealizable_ = true;
                subsetRealizable_ = true;
                nRealizableMoments_ = nMoments_;
            }
            else
            {
                checkCanonicalMoments(zeta, nN, nRealizableMoments_);

                if (nRealizableMoments_ == nMoments_)
                {
                    fullyRealizable_ = true;
                    subsetRealizable_ = true;
                }
                else if (nRealizableMoments_ > 1)
                {
                    fullyRealizable_ = false;
                    subsetRealizable_ = true;
                }
                // TODO Check if case of fully non-realizable is possible
            }

            calcNInvertibleMoments();
            realizabilityChecked_ = true;

            return;
        }
    }
}

void Foam::univariateMomentSet::calcNInvertibleMoments()
{
    if (nRealizableMoments_ % 2 != 0)
    {
        nInvertibleMoments_ = nRealizableMoments_ - 1;
    }
    else
    {
        nInvertibleMoments_ = nRealizableMoments_;
    }
}

void Foam::univariateMomentSet::setupQuadrature(bool clear)
{
    if (!realizabilityChecked_)
    {
        checkRealizability();
    }

    nNodes_ = nInvertibleMoments_/2.0;

    if (clear)
    {
        weights_.clear();
        abscissae_.clear();
    }

    weights_.resize(nNodes_, 0.0);
    abscissae_.resize(nNodes_, 0.0);
}

void Foam::univariateMomentSet::update()
{
    // NOTE Recomputing all the moments (even if they originally were
    //      not realizable) from quadrature.
    for (label momentI = 0; momentI < nMoments_; momentI++)
    {
        (*this)[momentI] = 0.0;

        for (label nodeI = 0; nodeI < nNodes_; nodeI++)
        {
            (*this)[momentI] += weights_[nodeI]*pow(abscissae_[nodeI], momentI);
        }
    }

    realizabilityChecked_ = false;
    inverted_ = false;
}

// ************************************************************************* //
