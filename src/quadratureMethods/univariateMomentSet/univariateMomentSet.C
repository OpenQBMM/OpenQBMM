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
    const word& quadratureType,
    const word& support,
    const scalar knownAbscissa
)
:
    scalarDiagonalMatrix(nMoments, initValue),
    nMoments_(nMoments),
    alpha_(),
    beta_(),
    quadratureType_(quadratureType),
    support_(support),
    negativeZeta_(0),
    degenerate_(false),
    inverted_(false),
    fullyRealizable_(true),
    subsetRealizable_(true),
    onMomentSpaceBoundary_(false),
    realizabilityChecked_(false),
    quadratureSetUp_(false),
    forceGauss_(false),
    nInvertibleMoments_(nMoments_),
    nRealizableMoments_(0),
    nNodes_(0),
    knownAbscissa_(knownAbscissa),
    weights_(),
    abscissae_()
{
    if (support_ != "R" && support_ != "RPlus" && support_ != "01")
    {
        FatalErrorInFunction
            << "The specified support is invalid." << nl
            << "    Valid supports are: R, RPlus and 01." << nl
            << "    Moment set: " << (*this)
            << abort(FatalError);
    }

    if (quadratureType_ != "Gauss" && quadratureType_ != "GaussRadau")
    {
        FatalErrorInFunction
            << "The specified quadrature type is invalid." << endl
            << "Valid quadrature types are: Gauss and GaussRadau."
            << abort(FatalError);
    }

    label recurrenceSize = label((nMoments_ - 2)/2) + 1;

    if (quadratureType_ == "GaussRadau")
    {
        recurrenceSize += 1;
    }

    alpha_.setSize(recurrenceSize, scalar(0));
    beta_.setSize(recurrenceSize + 1, scalar(0));
}

Foam::univariateMomentSet::univariateMomentSet
(
    const scalarDiagonalMatrix& m,
    const word& quadratureType,
    const word& support,
    const scalar knownAbscissa
)
:
    scalarDiagonalMatrix(m),
    nMoments_(m.size()),
    alpha_(),
    beta_(),
    quadratureType_(quadratureType),
    support_(support),
    negativeZeta_(0),
    degenerate_(false),
    inverted_(false),
    fullyRealizable_(true),
    subsetRealizable_(true),
    onMomentSpaceBoundary_(false),
    realizabilityChecked_(false),
    quadratureSetUp_(false),
    forceGauss_(false),
    nInvertibleMoments_(nMoments_),
    nRealizableMoments_(0),
    nNodes_(0),
    knownAbscissa_(knownAbscissa),
    weights_(),
    abscissae_()
{
    if (support_ != "R" && support_ != "RPlus" && support_ != "01")
    {
        FatalErrorInFunction
            << "The specified support is invalid." << nl
            << "    Valid supports are: R, RPlus and 01."
            << abort(FatalError);
    }

    if (quadratureType_ != "Gauss" && quadratureType_ != "GaussRadau")
    {
        FatalErrorInFunction
            << "The specified quadrature type is invalid." << endl
            << "Valid quadrature types are: Gauss and GaussRadau."
            << abort(FatalError);
    }

    label recurrenceSize = label((nMoments_ - 2)/2) + 1;

    alpha_.setSize(recurrenceSize, scalar(0));
    beta_.setSize(recurrenceSize + 1, scalar(0));
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

    if ((*this)[0] < SMALL)
    {
        nNodes_ = 0;

        return;
    }

    if (isDegenerate()) // || (*this)[0] <= 0.0)
    {
        //degenerate_ = true;
        setupQuadrature();
        weights_[0] = (*this)[0];
        abscissae_[0] = 0.0;
        inverted_ = true;

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
        FatalErrorInFunction
            << "Insufficient number (" << nInvertibleMoments_
            << ") of moments to define quadrature." << nl
            << "    Moment set: " << (*this)
            << abort(FatalError);
    }

    if (forceGauss_)
    {
        WarningInFunction
            << "Forcing Gauss quadrature. " << nl
            << "    Originally requested quadrature type: "
            << quadratureType_ << nl
            << "    Number of realizable moments: "
            << nRealizableMoments_ << nl
            << "    Moment set: " << (*this)
            << endl;
    }

    if (nInvertibleMoments_ == 2)
    {
        weights_[0] = (*this)[0];
        abscissae_[0] = (*this)[1]/(*this)[0];

        inverted_ = true;

        return;
    }

    if (quadratureType_ == "GaussRadau" && !forceGauss_)
    {
        // Compute P_{N-1} and P_{N-2} by recurrence
        // It is assumed the added point has abscissa in zero (xi0 = 0)
        scalar p = knownAbscissa_ - alpha_[0];
        scalar pMinus1 = 1.0;
        scalar p1 = p;

        for (label i = 1; i < nNodes_ - 1; i++)
        {
            p = knownAbscissa_ - alpha_[0]*p1 - beta_[i]*pMinus1;
            pMinus1 = p1;
            p1 = p;
        }

        alpha_[nNodes_ - 1] = knownAbscissa_ - beta_[nNodes_ - 1]*pMinus1/p;
    }

    scalarSquareMatrix z(nNodes_, scalar(0));

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
    const label nZeta
)
{
    scalarDiagonalMatrix canonicalMoments(nZeta, 0.0);

    canonicalMoments[0] = zeta[0];

    if (canonicalMoments[0] == 1.0)
    {
        nRealizableMoments_ = 2;
        onMomentSpaceBoundary_ = true;

        return;
    }

    for (label zetai = 1; zetai < nZeta; zetai++)
    {
        canonicalMoments[zetai] = zeta[zetai]/(1.0 - canonicalMoments[zetai-1]);

        if (canonicalMoments[zetai] < 0.0 || canonicalMoments[zetai] > 1.0)
        {
            nRealizableMoments_ = zetai + 1;

            return;
        }
        else if
        (
            canonicalMoments[zetai] == 0.0
         || canonicalMoments[zetai] == 1.0
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

void Foam::univariateMomentSet::checkRealizability()
{
    if (realizabilityChecked_)
    {
        return;
    }

    // If the zero-order moment is negative, exit immediately.
    if ((*this)[0] < 0.0)
    {
        FatalErrorInFunction
            << "The zero-order moment is negative." << nl
            << "    Moment set: " << (*this)
            << abort(FatalError);
    }

    // Check for the degenerate case where only m0 is defined
    if (nMoments_ <= 1)
    {
        FatalErrorInFunction
            << "The moment has size less or equal to 1." << nl
            << "    Moment set: " << (*this)
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

        if (zeta[0] <= 0.0)
        {
            if (isDegenerate() || zeta[0] == 0.0)
            {
                negativeZeta_ = 0;
                nRealizableMoments_ = 2;
                nInvertibleMoments_ = 2;
                fullyRealizable_ = true;
                subsetRealizable_ = true;
                onMomentSpaceBoundary_ = true;

                return;
            }

            negativeZeta_ = 1;
            nRealizableMoments_ = 1;
            nInvertibleMoments_ = 0;
            fullyRealizable_ = false;
            subsetRealizable_ = true;
            onMomentSpaceBoundary_ = false;

            FatalErrorInFunction
                << "Moment set with dimension 2 and only one realizable moment."
                << nl << "    Moment set: " << (*this)
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
            onMomentSpaceBoundary_ = false;

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

                if (zeta[0] < 1.0)
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
                    negativeZeta_ = 0;
                    nRealizableMoments_ = 2;
                    nInvertibleMoments_ = 2;
                    fullyRealizable_ = true;
                    subsetRealizable_ = true;
                    onMomentSpaceBoundary_ = true;

                    return;
                }

                negativeZeta_ = 1;
                nRealizableMoments_ = 1;
                nInvertibleMoments_ = 0;
                fullyRealizable_ = false;
                subsetRealizable_ = true;

                FatalErrorInFunction
                    << "Moment set with dimension 2 and only one "
                    << "realizable moment." << nl
                    << "    Moment set: " << (*this)
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

    if (!(support_ == "R") && zeta[0] <= 0.0)
    {
        if (isDegenerate() || zeta[0] == 0.0)
        {
            negativeZeta_ = 0;
            nRealizableMoments_ = 2;
            nInvertibleMoments_ = 2;
            fullyRealizable_ = false;
            subsetRealizable_ = true;
            onMomentSpaceBoundary_ = true;

            return;
        }

        negativeZeta_ = 1;
        nRealizableMoments_ = 1;
        nInvertibleMoments_ = 0;
        fullyRealizable_ = false;
        subsetRealizable_ = true;
        onMomentSpaceBoundary_ = false;

        FatalErrorInFunction
            << "Moment set with only one realizable moment." << nl
            << "    Moment set: " << (*this)
            << abort(FatalError);
    }

    for (label zetai = 1; zetai <= nD - 1; zetai++)
    {
        beta_[zetai] = zRecurrence[zetai][zetai]
                /zRecurrence[zetai - 1][zetai - 1];

        if (support_ == "R")
        {
            if (beta_[zetai] < 0.0)
            {
                nRealizableMoments_ = 2*zetai;
                calcNInvertibleMoments();
                fullyRealizable_ = false;
                subsetRealizable_ = true;
                realizabilityChecked_ = true;

                return;
            }
        }
        else
        {
            zeta[2*zetai - 1] = beta_[zetai]/zeta[2*zetai - 2];

            if (zeta[2*zetai - 1] <= 0.0)
            {
                if (support_ == "RPlus")
                {
                    if (zeta[2*zetai - 1] < 0.0)
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
                    checkCanonicalMoments(zeta, 2*zetai);
                }

                calcNInvertibleMoments();
                fullyRealizable_ = false;
                subsetRealizable_ = true;
                realizabilityChecked_ = true;

                return;
            }
        }

        alpha_[zetai] = zRecurrence[zetai][zetai + 1]/zRecurrence[zetai][zetai]
                - zRecurrence[zetai - 1][zetai]
                /zRecurrence[zetai - 1][zetai - 1];

        if (!(support_ == "R"))
        {
            zeta[2*zetai] = alpha_[zetai] - zeta[2*zetai - 1];

            if (zeta[2*zetai] <= 0.0)
            {
                if (support_ == "RPlus")
                {
                    if (zeta[2*zetai] < 0.0)
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
                    checkCanonicalMoments(zeta, 2*zetai + 1);
                }

                calcNInvertibleMoments();
                fullyRealizable_ = false;
                subsetRealizable_ = true;
                realizabilityChecked_ = true;

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

        if (zeta[2*nD - 1] <= 0.0)
        {
            if (support_ == "RPlus")
            {
                if (zeta[2*nD - 1] < 0.0)
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
                checkCanonicalMoments(zeta, 2*nD);
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

            if (zeta[2*nD] <= 0.0)
            {
                if (support_ == "RPlus")
                {
                    if (zeta[2*nD] < 0.0)
                    {
                        negativeZeta_ = 2*nD + 1;
                        nRealizableMoments_ = negativeZeta_;
                        fullyRealizable_ = false;
                        onMomentSpaceBoundary_ = false;
                    }
                    else
                    {
                        negativeZeta_ = 0;
                        nRealizableMoments_ = nMoments_;
                        fullyRealizable_ = true;
                        onMomentSpaceBoundary_ = true;
                    }
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta, 2*nD + 1);

                    if (onMomentSpaceBoundary_)
                    {
                        fullyRealizable_ = true;
                    }
                    else
                    {
                        fullyRealizable_ = false;
                    }
                }

                calcNInvertibleMoments();

                subsetRealizable_ = true;
                realizabilityChecked_ = true;

                return;
            }
            else // zeta[2*nD] > 0.0
            {
                if (support_ == "RPlus")
                {
                    negativeZeta_ = 0;
                    nRealizableMoments_ = nMoments_;
                    fullyRealizable_ = true;
                    onMomentSpaceBoundary_ = false;
                }
                else // Support on [0,1]
                {
                    checkCanonicalMoments(zeta, 2*nD + 1);

                    if (onMomentSpaceBoundary_)
                    {
                        fullyRealizable_ = true;
                    }
                    else
                    {
                        fullyRealizable_ = false;
                    }
                }

                calcNInvertibleMoments();

                subsetRealizable_ = true;
                realizabilityChecked_ = true;

                return;
            }
        }
        else
        {
            // If support is [0, + inf[ and this level is reached, the full set
            // of moments is realizable
            if (support_ == "RPlus")
            {
                fullyRealizable_ = true;
                subsetRealizable_ = true;
                nRealizableMoments_ = nMoments_;
                onMomentSpaceBoundary_ = false;
            }
            else
            {
                checkCanonicalMoments(zeta, nN);

                if (nRealizableMoments_ == nMoments_)
                {
                    fullyRealizable_ = true;
                    subsetRealizable_ = true;
                }
                else
                {
                    fullyRealizable_ = false;
                    subsetRealizable_ = true;
                }
            }

            calcNInvertibleMoments();
            realizabilityChecked_ = true;

            return;
        }
    }
}

void Foam::univariateMomentSet::calcNInvertibleMoments()
{
    if (quadratureType_ == "Gauss")
    {
        if (nRealizableMoments_ % 2 != 0)
        {
            nInvertibleMoments_ = nRealizableMoments_ - 1;
        }
        else
        {
            nInvertibleMoments_ = nRealizableMoments_;
        }

        return;
    }
    else
    {
        forceGauss_ = false;
        nInvertibleMoments_ = nRealizableMoments_;

        if (nRealizableMoments_ % 2 == 0)
        {
            forceGauss_ = true;
        }
    }
}

void Foam::univariateMomentSet::setupQuadrature
(
    bool clear,
    bool nullMomentSet
)
{
    if (!realizabilityChecked_)
    {
        checkRealizability();
    }

    if (degenerate_)
    {
        nNodes_ = 1;
    }
    else
    {
        nNodes_ = nInvertibleMoments_/2;

        if (quadratureType_ == "GaussRadau" && !forceGauss_)
        {
            nNodes_ += 1;
        }
    }

    if (clear)
    {
        weights_.clear();
        abscissae_.clear();
    }

    weights_.resize(nNodes_, 0.0);
    abscissae_.resize(nNodes_, 0.0);

    //quadratureSetUp_ = true;
}

void Foam::univariateMomentSet::update()
{
    // Recomputing all the moments (even if they originally were not realizable)
    // from quadrature (projection step).
    for (label momenti = 0; momenti < nMoments_; momenti++)
    {
        (*this)[momenti] = 0.0;

        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            (*this)[momenti] += weights_[nodei]*pow(abscissae_[nodei], momenti);
        }
    }

    realizabilityChecked_ = false;
    inverted_ = false;
}

// ************************************************************************* //
