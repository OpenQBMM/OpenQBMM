/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 Alberto Passalacqua
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
    const scalar initValue
)
:
    scalarDiagonalMatrix(nMoments, initValue),
    nMoments_(nMoments),
    inverted_(false),
    realizable_(true),
    realizabilityChecked_(false),
    quadratureSetUp_(false),
    nInvertibleMoments_(nMoments_)
{}

Foam::univariateMomentSet::univariateMomentSet
(
    const scalarDiagonalMatrix& m
)
:
    scalarDiagonalMatrix(m),
    nMoments_(m.size()),
    inverted_(false),
    realizable_(true),
    realizabilityChecked_(false),
    quadratureSetUp_(false),
    nInvertibleMoments_(nMoments_)
{}


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
        isRealizable();
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

    // Storage for the recurrence relationship
    scalarDiagonalMatrix a(nNodes_, scalar(0));
    scalarDiagonalMatrix b(nNodes_, scalar(0));
    scalarSquareMatrix sig(2*nNodes_ + 1, 2*nNodes_ + 1, scalar(0));

    // Applying Wheeler's algorithm
    for (label i = 1; i < nInvertibleMoments_ + 1; i++)
    {
        sig[1][i] = (*this)[i-1];
    }

    a[0] = (*this)[1]/(*this)[0];
    b[0] = scalar(0);

    for (label k = 2; k < nNodes_ + 1; k++)
    {
        for (label l = k; l < nInvertibleMoments_ - k + 3; l++)
        {
            sig[k][l] = sig[k-1][l+1] - a[k-2]*sig[k-1][l] - b[k-2]*sig[k-2][l];
        }

        a[k-1] = sig[k][k+1]/sig[k][k] - sig[k-1][k]/sig[k-1][k-1];
        b[k-1] = sig[k][k]/sig[k-1][k-1] ;
    }

    // Checking moment realizability
    for (label i = 0; i < nNodes_; i++)
    {
        if (b[i] < scalar(0))
        {
            FatalErrorIn
            (
                "Foam::univariateMomentSet::invert\n"
                "(\n"
                "    const scalarDiagonalMatrix& weights,\n"
                "    const scalarDiagonalMatrix& abscissae\n"
                ")"
            )   << "Moments are not realizable."
                << abort(FatalError);
        }
    }

    scalarSquareMatrix z(nNodes_, nNodes_, scalar(0));

    for (label i = 0; i < nNodes_ - 1; i++)
    {
        z[i][i] = a[i];
        z[i][i+1] = Foam::sqrt(b[i+1]);
        z[i+1][i] = z[i][i+1];
    }

    z[nNodes_ - 1][nNodes_ - 1] = a[nNodes_ - 1];

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

void Foam::univariateMomentSet::checkRealizability() 
{
    if (realizabilityChecked_)
    {
        return;
    }
    
    label nN = nMoments_ - 1;
    label nD = label(nN/2);
    label nR = nN - 2*nD;

    label nZRows = nD + 1;
    label nZColumns = nMoments_;
    label nAlpha = label((nN - 1)/2) + 1;
    label nBeta = nD + 1;

    scalarRectangularMatrix z(nZRows, nZColumns, 0.0);
    scalarDiagonalMatrix alpha(nAlpha, 0.0);
    scalarDiagonalMatrix beta(nBeta, 0.0);
    scalarDiagonalMatrix zeta(nN);

    // Check for the degenerate case where only m0 is defined
    if (nMoments_ <= 1)
    {
        FatalErrorIn
        (
            "Foam::univariateMomentSet::isRealizable()\n"
        )   << "The moment set is degenerate."
            << abort(FatalError);
    }

    // Check if m0 < 0
    if ((*this)[0] < 0)
    {
        FatalErrorIn
        (
            "Foam::univariateMomentSet::isRealizable()\n"
        )   << "The zero-order moment is negative."
            << abort(FatalError);
    }

    // Check for the case with only two moments
    if (nMoments_ == 2)
    {
        zeta[0] = (*this)[1]/(*this)[0] - 1.0;

        if (zeta[0] <= 0)
        {
            negativeZeta_ = 1;
            nRealizableMoments_ = 1;
            nInvertibleMoments_ = 0;
            realizable_ = false;

            FatalErrorIn
            (
                "Foam::univariateMomentSet::isRealizable()\n"
            )   << "Moment set with dimension 2 and only one realizable moment."
                << abort(FatalError);
        }

        negativeZeta_ = 0;
        nRealizableMoments_ = 2;
        nInvertibleMoments_ = 2;
        realizable_ = true;

        return;
    }

    // Fill the first row of the z matrix with the moments
    for (label columnI = 0; columnI < nZColumns; columnI++)
    {
        z[0][columnI] = (*this)[columnI];
    }

    alpha[0] = (*this)[1]/(*this)[0];
    beta[0] = (*this)[0];

    for (label columnI = 1; columnI < nZColumns - 1; columnI++)
    {
        z[1][columnI] = z[0][columnI + 1] - alpha[0]*z[0][columnI];
    }

    zeta[0] = alpha[0];

    if (zeta[0] <= 0)
    {
        negativeZeta_ = 1;
        nRealizableMoments_ = 1;
        nInvertibleMoments_ = 0;
        realizable_ = false;
        
        FatalErrorIn
        (
            "Foam::univariateMomentSet::isRealizable()\n"
        )   << "Moment set with only one realizable moment."
            << abort(FatalError);
    }

    for (label zetaI = 1; zetaI <= nD - 1; zetaI++)
    {
        beta[zetaI] = z[zetaI][zetaI]/z[zetaI - 1][zetaI - 1];
        zeta[2*zetaI - 1] = beta[zetaI]/zeta[2*zetaI - 2];

        if (zeta[2*zetaI - 1] <= 0)
        {
            negativeZeta_ = 2*zetaI;
            nRealizableMoments_ = negativeZeta_;
            nInvertibleMoments_ = nRealizableMoments_;
            realizable_ = false;

            return;
        }

        alpha[zetaI] = z[zetaI][zetaI + 1]/z[zetaI][zetaI] 
                - z[zetaI - 1][zetaI]/z[zetaI - 1][zetaI - 1];

        zeta[2*zetaI] = alpha[zetaI] - zeta[2*zetaI - 1];

        if (zeta[2*zetaI] <= 0)
        {
            negativeZeta_ = 2*zetaI + 1;
            nRealizableMoments_ = negativeZeta_;
            nInvertibleMoments_ = nRealizableMoments_ - 1;
            realizable_ = false;
    
            return;
        }

        for (label columnI = zetaI + 1; columnI <= nN - zetaI - 1; columnI++)
        {
            z[zetaI + 1][columnI] = z[zetaI][columnI + 1] 
                    - alpha[zetaI]*z[zetaI][columnI]
                    - beta[zetaI]*z[zetaI - 1][columnI];
        }
    }

    beta[nD] = z[nD][nD]/z[nD - 1][nD - 1];
    zeta[2*nD - 1] = beta[nD]/zeta[2*nD - 2];

    if (zeta[2*nD - 1] <= 0)
    {
        negativeZeta_ = 2*nD;
        nRealizableMoments_ = negativeZeta_;
        nInvertibleMoments_ = nRealizableMoments_;
        realizable_ = false;

        return;
    }

    if (nR == 1)
    {
        alpha[nD] = z[nD][nD + 1]/z[nD][nD] - z[nD - 1][nD]/z[nD - 1][nD - 1];
        zeta[2*nD] = alpha[nD] - zeta[2*nD - 1];

        if (zeta[2*nD] <= 0)
        {
            negativeZeta_ = 2*nD + 1;
            nRealizableMoments_ = negativeZeta_;
            nInvertibleMoments_ = nRealizableMoments_ - 1;
            realizable_ = false;

            return;
        }
    }

    realizable_ = true;
    nRealizableMoments_ = nMoments_;

    if (nMoments_ % 2 != 0)
    {
        nInvertibleMoments_ = nMoments_ - 1;
    }
    else
    {
        nInvertibleMoments_ = nMoments_;
    }
}

void Foam::univariateMomentSet::setupQuadrature(bool clear)
{
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
    //      Should we limit to nRealizableMoments_?
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
