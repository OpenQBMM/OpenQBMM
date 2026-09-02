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

#include "gaussLobattoMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gaussLobattoMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentInversion,
        gaussLobattoMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gaussLobattoMomentInversion::gaussLobattoMomentInversion
(
    const dictionary& dict,
    const label nMaxNodes
)
:
    univariateMomentInversion(dict, nMaxNodes),
    forceRadau_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gaussLobattoMomentInversion::correctRecurrence
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    if (forceRadau_)
    {
        // Only the lower abscissa is fixed, as in Gauss-Radau
        scalar pMinus1 = 1.0;

        const scalar p =
            orthogonalPolynomial(alpha, beta, minKnownAbscissa, pMinus1);

        alpha[nNodes_ - 1] =
                minKnownAbscissa - beta[nNodes_ - 1]*pMinus1/p;
    }
    else
    {
        // Both abscissae are fixed, which requires correcting the last alpha
        // and beta coefficients by solving a 2x2 system
        scalar pMinus1Left = 1.0;
        scalar pMinus1Right = 1.0;

        const scalar pLeft =
            orthogonalPolynomial(alpha, beta, minKnownAbscissa, pMinus1Left);

        const scalar pRight =
            orthogonalPolynomial(alpha, beta, maxKnownAbscissa, pMinus1Right);

        const scalar d = pLeft*pMinus1Right - pRight*pMinus1Left;

        alpha[nNodes_ - 1] =
                (minKnownAbscissa*pLeft*pMinus1Right
                - maxKnownAbscissa*pRight*pMinus1Left)/d;

        beta[nNodes_ - 1] =
                (maxKnownAbscissa - minKnownAbscissa)*pLeft*pRight/d;
    }
}

void Foam::gaussLobattoMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    label nRealizableMoments = moments.nRealizableMoments();

    if (nRealizableMoments > 2)
    {
        if (nRealizableMoments % 2 == 0)
        {
            nInvertibleMoments_ = nRealizableMoments;
            forceRadau_ = false;
            nNodes_ = nInvertibleMoments_/2 + 1;
        }
        else
        {
            nInvertibleMoments_ = nRealizableMoments;
            forceRadau_ = true;
            nNodes_ = nInvertibleMoments_/2 + 1;
        }
    }
    else
    {
        FatalErrorInFunction
            << "The moment set has two or less realizable moments." << nl
            << "    Moment set: " << moments
            << abort(FatalError);
    }

    abscissae_.setSize(nNodes_);
    weights_.setSize(nNodes_);
}


Foam::label Foam::gaussLobattoMomentInversion::nAdditionalQuadraturePoints
(
    const label nMoments
) const
{
    return 2;
}

// ************************************************************************* //
