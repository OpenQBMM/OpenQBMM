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
    realizable_(true),
    realizabilityChecked_(false),
    nInvertibleMoments_(nMoments_)
{};

Foam::univariateMomentSet::univariateMomentSet
(
    const scalarDiagonalMatrix& m
)
:
    scalarDiagonalMatrix(m),
    nMoments_(m.size()),
    realizable_(true),
    realizabilityChecked_(false),
    nInvertibleMoments_(nMoments_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateMomentSet::~univariateMomentSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::univariateMomentSet::invert
(
    scalarDiagonalMatrix& weights,
    scalarDiagonalMatrix& abscissae
)
{
    if (!realizabilityChecked_)
    {
        isRealizable();
    }

    if (nInvertibleMoments_ < 2)
    {
        //Info << "Moments: " << (*this) << endl;
        
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

    // Number of quadrature nodes based on nInvertibleMoments_ (always even)
    label nNodes = nInvertibleMoments_/2.0;

    // Storage for the recurrence relationship
    scalarDiagonalMatrix a(nNodes, scalar(0));
    scalarDiagonalMatrix b(nNodes, scalar(0));
    scalarSquareMatrix sig(2*nNodes + 1, 2*nNodes + 1, scalar(0));

    // Applying Wheeler's algorithm
    for (label i = 1; i < nInvertibleMoments_ + 1; i++)
    {
        sig[1][i] = (*this)[i-1];
    }

    a[0] = (*this)[1]/(*this)[0];
    b[0] = scalar(0);

    for (label k = 2; k < nNodes + 1; k++)
    {
        for (label l = k; l < nInvertibleMoments_ - k + 3; l++)
        {
            sig[k][l] = sig[k-1][l+1] - a[k-2]*sig[k-1][l] - b[k-2]*sig[k-2][l];
        }

        a[k-1] = sig[k][k+1]/sig[k][k] - sig[k-1][k]/sig[k-1][k-1];
        b[k-1] = sig[k][k]/sig[k-1][k-1] ;
    }

    // Checking moment realizability
    for (label i = 0; i < nNodes; i++)
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

    scalarSquareMatrix z(nNodes, nNodes, scalar(0));

    for (label i = 0; i < nNodes - 1; i++)
    {
        z[i][i] = a[i];
        z[i][i+1] = Foam::sqrt(b[i+1]);
        z[i+1][i] = z[i][i+1];
    }

    z[nNodes - 1][nNodes - 1] = a[nNodes - 1];

    // Computing weights and abscissae
    eigenSolver zEig(z, true);

    // Resetting previous weights and abscissae to zero. Needed due to 
    // adaption performed through the realizability check
    forAll (weights, nodeI)
    {
        weights[nodeI] = 0.0;
        abscissae[nodeI] = 0.0;
    }

    // Storing weight and abscissae in destination containers.
    // Only computed values are copied. The remaining ones are left equal to
    // zero.
    for (label i = 0; i < nNodes; i++)
    {
        weights[i] = (*this)[0]*sqr(zEig.eigenvectors()[0][i]);
        abscissae[i] = zEig.eigenvaluesRe()[i];
    }
}

bool Foam::univariateMomentSet::isRealizable() 
{
    if (realizabilityChecked_)
    {
    return realizable_;
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
            nInvertibleMoments_ = 0;
            realizable_ = false;

            FatalErrorIn
            (
                "Foam::univariateMomentSet::isRealizable()\n"
            )   << "Moment set with dimension 2 and only one realizable moment."
                << abort(FatalError);
        }

        negativeZeta_ = 0;
        nInvertibleMoments_ = 2;
        realizable_ = true;

        return realizable_;
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
        nInvertibleMoments_ = 1;
        realizable_ = false;

        return realizable_;
    }

    for (label zetaI = 1; zetaI <= nD - 1; zetaI++)
    {
        beta[zetaI] = z[zetaI][zetaI]/z[zetaI - 1][zetaI - 1];
        zeta[2*zetaI - 1] = beta[zetaI]/zeta[2*zetaI - 2];

        if (zeta[2*zetaI - 1] <= 0)
        {
            negativeZeta_ = 2*zetaI;
            nInvertibleMoments_ = negativeZeta_;
            realizable_ = false;

            return realizable_;
        }

        alpha[zetaI] = z[zetaI][zetaI + 1]/z[zetaI][zetaI] 
                - z[zetaI - 1][zetaI]/z[zetaI - 1][zetaI - 1];

        zeta[2*zetaI] = alpha[zetaI] - zeta[2*zetaI - 1];

        if (zeta[2*zetaI] <= 0)
        {
            negativeZeta_ = 2*zetaI + 1;  
            nInvertibleMoments_ = negativeZeta_ - 1;
            realizable_ = false;
    
            return realizable_;
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
        nInvertibleMoments_ = negativeZeta_;
        realizable_ = false;

        return realizable_;
    }

    if (nR == 1)
    {
        alpha[nD] = z[nD][nD + 1]/z[nD][nD] - z[nD - 1][nD]/z[nD - 1][nD - 1];
        zeta[2*nD] = alpha[nD] - zeta[2*nD - 1];

        if (zeta[2*nD] <= 0)
        {
            negativeZeta_ = 2*nD + 1;
            nInvertibleMoments_ = negativeZeta_ - 1;
            realizable_ = false;

            return realizable_;
        }
    }

    realizable_ = true;

    if (nMoments_ % 2 != 0)
    {
        nInvertibleMoments_ = nMoments_ - 1;
    }

    return realizable_;
}

void Foam::univariateMomentSet::update
(
    const scalarDiagonalMatrix& weights, 
    const scalarDiagonalMatrix& abscissae
)
{
    for (label momentI = 0; momentI < nMoments_; momentI++)
    {
    (*this)[momentI] = 0.0;

    for (label nodeI = 0; nodeI < weights.size(); nodeI++)
    {
        (*this)[momentI] += weights[nodeI]*pow(abscissae[nodeI], momentI);
    }
    }

    realizabilityChecked_ = false;
}

// ************************************************************************* //
