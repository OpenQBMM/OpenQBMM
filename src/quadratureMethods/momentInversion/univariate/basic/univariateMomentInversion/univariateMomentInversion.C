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
    Copyright (C) 2019-2020 Alberto Passalacqua
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
#include "EigenMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(univariateMomentInversion, 0);
    defineRunTimeSelectionTable(univariateMomentInversion, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMomentInversion::univariateMomentInversion
(
    const dictionary& dict
)
:
    smallM0_(dict.lookupOrDefault<scalar>("smallM0", 1.0e-12)),
    nInvertibleMoments_(),
    nNodes_(),
    abscissae_(),
    weights_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateMomentInversion::~univariateMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::univariateMomentInversion::JacobiMatrix
(
    univariateMomentSet& moments,
    scalarSquareMatrix& z,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    scalarList alpha(moments.alphaRecurrence());
    scalarList beta(moments.betaRecurrence());

    correctRecurrence
    (
        moments,
        alpha,
        beta,
        minKnownAbscissa,
        maxKnownAbscissa
    );

    for (label i = 0; i < nNodes_ - 1; i++)
    {
        z[i][i] = alpha[i];
        z[i][i+1] = Foam::sqrt(beta[i+1]);
        z[i+1][i] = z[i][i+1];
    }

    z[nNodes_ - 1][nNodes_ - 1] = alpha[nNodes_ - 1];
}

void Foam::univariateMomentInversion::invert
(
    univariateMomentSet& moments,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    if (moments.isDegenerate())
    {
        nNodes_ = 1;
        weights_.setSize(nNodes_);
        abscissae_.setSize(nNodes_);
        weights_[0] = moments[0];
        abscissae_[0] = 0.0;

        return;
    }

    if (moments[0] < smallM0_)
    {
        nNodes_ = 0;

        weights_.setSize(nNodes_);
        abscissae_.setSize(nNodes_);

        return;
    }

    calcNQuadratureNodes(moments);

    if (nInvertibleMoments_ == 2)
    {
        weights_[0] = moments[0];
        abscissae_[0] = moments[1]/moments[0];

        return;
    }

    scalarSquareMatrix z(nNodes_, scalar(0));
    JacobiMatrix(moments, z, minKnownAbscissa, maxKnownAbscissa);

    // Computing weights and abscissae
    EigenMatrix<scalar> zEig(z, true);

    // Computing weights and abscissae
    for (label i = 0; i < nNodes_; i++)
    {
        weights_[i] = moments[0]*sqr(zEig.EVecs()[0][i]);
        abscissae_[i] = zEig.EValsRe()[i];
    }
}


// ************************************************************************* //
