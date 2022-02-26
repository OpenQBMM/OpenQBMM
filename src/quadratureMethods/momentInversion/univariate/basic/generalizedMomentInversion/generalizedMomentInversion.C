/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 Alberto Passalacqua
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

#include "generalizedMomentInversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(generalizedMomentInversion, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentInversion,
        generalizedMomentInversion,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::generalizedMomentInversion::generalizedMomentInversion
(
    const dictionary& dict,
    const label nMaxNodes
)
:
    univariateMomentInversion(dict, nMaxNodes),
    a_
    (
        dict.lookupOrDefault<scalar>
        (
            "a",
            0.0
        )
    ),
    nu_
    (
        dict.lookupOrDefault<scalar>
        (
            "nu",
            1.0
        )
    ),
    nMaxNodes_(nMaxNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::generalizedMomentInversion::~generalizedMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::generalizedMomentInversion::correctRecurrence
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    const word& support = moments.support();

    if (support == "R")
    {
        correctRecurrenceR(alpha, beta);
        return;
    }

    if (support == "RPlus")
    {
        correctRecurrenceRPlus(moments, alpha, beta);
        return;
    }

    if (support == "01")
    {
        correctRecurrence01(moments, alpha, beta);
        return;
    }
}

void Foam::generalizedMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    label nRealizableMoments = moments.nRealizableMoments(); // Trigger calculations of zeta_k
    
    Info << "nRealizableMoments = " << nRealizableMoments << endl;
    nRegularQuadratureNodes_ 
        = (nRealizableMoments % 2 != 0) 
        ? label((nRealizableMoments - 1)/2.0) 
        : label(nRealizableMoments/2.0);

    if (nRealizableMoments > 3)
    {
        nAdditionalQuadratureNodes_ = nMaxNodes_ - nRegularQuadratureNodes_;
        nNodes_ = nMaxNodes_;
    }
    else
    {
        nAdditionalQuadratureNodes_ = 0;
        nNodes_ = nRegularQuadratureNodes_;
    }
    
    abscissae_.setSize(nMaxNodes_);
    weights_.setSize(nMaxNodes_);

    Info << "nNaxNodes = " << nMaxNodes_ << endl;
    Info << "nRegularQuadratureNodes = " 
         << nRegularQuadratureNodes_ << endl;
    Info << "nAdditionalQuadratureNodes = " 
         << nAdditionalQuadratureNodes_ << endl;
}

void Foam::generalizedMomentInversion::invert
(
    univariateMomentSet& moments,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    (*this).univariateMomentInversion::invert
        (
            moments, 
            minKnownAbscissa, 
            maxKnownAbscissa
        );
}

void Foam::generalizedMomentInversion::correctRecurrenceR
(
    scalarList& alpha,
    scalarList& beta
)
{
    scalar an = 0;
    
    for (label i = 0; i < nRegularQuadratureNodes_; i++)
    {
        an += alpha[i];
    }

    an /= nRegularQuadratureNodes_;

    for (label i = nRegularQuadratureNodes_; i < nNodes_; i++)
    {
        alpha[i] = an;

        beta[i-1] = beta[nRegularQuadratureNodes_ - 1]
                   *pow(scalar(i - 1)/scalar(nRegularQuadratureNodes_ -1), nu_);
    }

    beta[nNodes_ - 1] = beta[nRegularQuadratureNodes_ - 1]
                        *pow(scalar(nNodes_ - 1)
                        /scalar(nRegularQuadratureNodes_ -1), nu_);
}

void Foam::generalizedMomentInversion::correctRecurrenceRPlus
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta
)
{
    label nMoments = moments.size();
    scalarList& zetas(moments.zetas());

    zetas.resize(2*nMaxNodes_);

    if (nMoments + 1 <= 2*nMaxNodes_)
    {
        zetas[nMoments] = (nMoments + 1.0)*zetas[2*nRegularQuadratureNodes_ - 1]
                          /(2.0*nRegularQuadratureNodes_);
    }

    for 
    (
        label i = nRegularQuadratureNodes_; 
        i < nMaxNodes_ && nAdditionalQuadratureNodes_ > 0; 
        i++
    )
    {
        zetas[2*i - 1] = (i + a_)*zetas[2*nRegularQuadratureNodes_ - 3]
                        /(nRegularQuadratureNodes_ + a_ - 1.0);

        zetas[2*i] = (i + 1)*zetas[2*nRegularQuadratureNodes_ - 2]
                    /(nRegularQuadratureNodes_);
    }

    alpha[0] = zetas[0];

    for (label i = 1; i < nMaxNodes_; i++)
    {
        alpha[i] = zetas[2*i] + zetas[2*i - 1];
    }

    for (label i = 1; i < nMaxNodes_; i++)
    {
        beta[i] = zetas[2*i - 1]*zetas[2*i - 2];
    }
}

void Foam::generalizedMomentInversion::correctRecurrence01
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta
)
{
    scalarList& zetas(moments.zetas());
    scalarList& canonicalMoments(moments.canonicalMoments());

    zetas.resize(2*nMaxNodes_);
    canonicalMoments.resize(2*nMaxNodes_);

    for
    (
        label i = nRegularQuadratureNodes_; 
        i < nMaxNodes_ && nAdditionalQuadratureNodes_ > 0; 
        i++
    )
    {
        canonicalMoments[2*i - 1] = canonicalMoments[2*i - 3];
        canonicalMoments[2*i] = canonicalMoments[2*i - 2];

        zetas[2*i - 1] 
            = canonicalMoments[2*i - 1]*(1.0 - canonicalMoments[2*i - 2]);

        zetas[2*i] 
            = canonicalMoments[2*i]*(1.0 - canonicalMoments[2*i - 1]);
    }

    alpha[0] = zetas[0];

    for (label i = 1; i < nMaxNodes_; i++)
    {
        alpha[i] = zetas[2*i] + zetas[2*i - 1];
    }

    for (label i = 1; i < nMaxNodes_; i++)
    {
        beta[i] = zetas[2*i - 1]*zetas[2*i - 2];
    }
}

// ************************************************************************* //
