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

#include "univariateMomentInversion.H"

#include "eigenSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMomentInversion::univariateMomentInversion
(
    univariateMomentSet& moments
)
:
    moments_(moments),
    inverted_(false),
    nInvertibleMoments_(moments.size()),
    nNodes_(),
    weights_(),
    abscissae_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateMomentInversion::~univariateMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarList& Foam::univariateMomentInversion::alphaRecurrence()
{
    return moments_.alphaRecurrence();
}

Foam::scalarList& Foam::univariateMomentInversion::betaRecurrence()
{
    return moments_.betaRecurrence();
}

Foam::scalarSquareMatrix Foam::univariateMomentInversion::JacobiMatrix
(
    const scalarList& alpha,
    const scalarList& beta
)
{
    scalarSquareMatrix z(nNodes_, scalar(0));

    for (label i = 0; i < nNodes_ - 1; i++)
    {
        z[i][i] = moments_.alphaRecurrence()[i];
        z[i][i+1] = Foam::sqrt(moments_.betaRecurrence()[i+1]);
        z[i+1][i] = z[i][i+1];
    }

    z[nNodes_ - 1][nNodes_ - 1] = moments_.alphaRecurrence()[nNodes_ - 1];

    return z;
}

void Foam::univariateMomentInversion::invert()
{
    if (inverted_ || moments_.isDegenerate())
    {
        return;
    }

    if (moments_[0] < SMALL)
    {
        nNodes_ = 0;

        return;
    }

    if (nInvertibleMoments_ == 2)
    {
        weights_[0] = moments_[0];
        abscissae_[0] = moments_[1]/moments_[0];

        inverted_ = true;

        return;
    }

    scalarList& a(alphaRecurrence());
    scalarList& b(betaRecurrence());

    scalarSquareMatrix z(JacobiMatrix(a, b));

    // Computing weights and abscissae
    eigenSolver zEig(z, true);

    // Computing weights and abscissae
    for (label i = 0; i < nNodes_; i++)
    {
        weights_[i] = moments_[0]*sqr(zEig.eigenvectors()[0][i]);
        abscissae_[i] = zEig.eigenvaluesRe()[i];
    }

    inverted_ = true;
}


// ************************************************************************* //
