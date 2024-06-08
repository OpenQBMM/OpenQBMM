/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 Alberto Passalacqua
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
    nu_
    (
        dict.lookupOrDefault<scalar>("nu", 1.0)
    ),
    ndfTypeRPlus_
    (
        dict.lookupOrDefault<word>("ndfTypeRPlus", "gamma")
    ),
    nMaxNodes_(nMaxNodes)
{
    
    if ((ndfTypeRPlus_ != "gamma" && ndfTypeRPlus_ != "lognormal"))
    {
        FatalErrorInFunction
            << "The specified type of NDF for RPlus must be gamma or" << nl
            << "lognormal. The current value is " << ndfTypeRPlus_
            << abort(FatalError);
    }
}


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
    }
    else if (support == "RPlus")
    {
        correctRecurrenceRPlus(moments, alpha, beta);
    }
    else if (support == "01")
    {
        correctRecurrence01(moments, alpha, beta);
    }
}

void Foam::generalizedMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    // Trigger calculations of zeta_k by computing the number of realizable
    // moments.
    label nRealizableMoments = moments.nRealizableMoments(); 
    
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

    /*Info << "nNaxNodes = " << nMaxNodes_ << endl
         << "nRegularQuadratureNodes = " 
         << nRegularQuadratureNodes_ << endl
         << "nAdditionalQuadratureNodes = " 
         << nAdditionalQuadratureNodes_ << endl;*/
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
    // The realizability test will establish if adding nodes with GQMOM is
    // feasible and set nAdditionalQuadratureNodes_ if not.
    if (nAdditionalQuadratureNodes_ <= 0)
    {
        return; // Use Gauss if no additional nodes are possible 
    }

    scalar an = 0;
    
    for (label i = 0; i < nRegularQuadratureNodes_; i++)
    {
        an += alpha[i];
    }

    an /= nRegularQuadratureNodes_;

    for (label i = nRegularQuadratureNodes_; i < nNodes_; i++)
    {
        alpha[i] = an;

        beta[i-1] = beta[nRegularQuadratureNodes_ - 1]*pow(scalar(i - 1)
                   /scalar(nRegularQuadratureNodes_ - 1), nu_);
    }

    beta[nNodes_ - 1] = beta[nRegularQuadratureNodes_ - 1]
                       *pow(scalar(nNodes_ - 1)
                       /scalar(nRegularQuadratureNodes_ - 1), nu_);
}

void Foam::generalizedMomentInversion::correctRecurrenceRPlus
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta
)
{
    // The realizability test will establish if adding nodes with GQMOM is
    // feasible and set nAdditionalQuadratureNodes_ if not.
    if (nAdditionalQuadratureNodes_ <= 0)
    {
        return; // Use Gauss if no additional nodes are possible 
    }

    moments.zetas().resize(2*nMaxNodes_ - 1, 0.0);

    if (ndfTypeRPlus_ == "gamma")
    {
        const scalar m1sqr = sqr(moments(1));
        const scalar alphaCoeff = m1sqr/(moments(0)*moments(2) - m1sqr) - 1.0; 

        for 
        (
            label i = nRegularQuadratureNodes_; 
            i < nMaxNodes_ && nAdditionalQuadratureNodes_ > 0; 
            i++
        )
        {
            moments.zetas()[2*i - 1] = 
                (i + alphaCoeff)*moments.zetas()[2*nRegularQuadratureNodes_ - 3]
               /(nRegularQuadratureNodes_ - 1 + alphaCoeff);

            moments.zetas()[2*i] = 
                (i + 1)*moments.zetas()[2*nRegularQuadratureNodes_ - 2]
               /(nRegularQuadratureNodes_);

            // Debug info
            //Info << "zetas[2*i-1] = " << zetas[2*i - 1] << endl;
            //Info << "2i-1 = " << 2*i - 1 << endl;
            //Info << "zetas[2*i] = " << zetas[2*i] << endl;
            //Info << "2i = " << 2*i << endl;
        }   
    }
    else if (ndfTypeRPlus_ == "lognormal")
    {
        const scalar eta = sqrt(moments(0)*moments(2)/sqr(moments(1)));

        for
        (
            label i = nRegularQuadratureNodes_; 
            i < nMaxNodes_ && nAdditionalQuadratureNodes_ > 0; 
            i++
        )
        {
            moments.zetas()[2*i - 1] = 
                pow(eta, 2*(i + 1 - nRegularQuadratureNodes_))
               *(
                    (pow(eta, 2*(i+1)) - 1.0)
                   /(pow(eta, 2*nRegularQuadratureNodes_) - 1.0)
                )
               *moments.zetas()[2*nRegularQuadratureNodes_ - 3];

            moments.zetas()[2*i] = 
                pow(eta, 4*(i + 1 - nRegularQuadratureNodes_))
               *moments.zetas()[2*nRegularQuadratureNodes_ - 2];
            // Info << "zetas[2*i-1] = " << moments.zetas()[2*i - 1] << endl;
            // Info << "2i-1 = " << 2*i - 1 << endl;
            // Info << "zetas[2*i] = " << moments.zetas()[2*i] << endl;
            // Info << "2i = " << 2*i << endl;
        }
    }    

    alpha[0] = moments.zetas()[0];

    for (label i = 1; i < nMaxNodes_; i++)
    {
        alpha[i] = moments.zetas()[2*i] + moments.zetas()[2*i - 1];
    }

    for (label i = 1; i < nMaxNodes_; i++)
    {
        beta[i] = moments.zetas()[2*i - 1]*moments.zetas()[2*i - 2];
    }
}

void Foam::generalizedMomentInversion::correctRecurrence01
(
    univariateMomentSet& moments,
    scalarList& alpha,
    scalarList& beta
)
{
    // The realizability test will establish if adding nodes with GQMOM is
    // feasible and set nAdditionalQuadratureNodes_ = 0 if not.
    if (nAdditionalQuadratureNodes_ <= 0)
    {
        return; // Use Gauss if no additional nodes are possible 
    }
    
    // We do not store z0 = 1, so we have 2*nRegularQuadratureNodes_ - 1 zetas
    moments.zetas().resize(2*nMaxNodes_ - 1);

    // We do not store p0, so canonicalMoments[0] = p1, which means that
    // we have 2*nRegularQuadratureNodes_ - 1 canonical moments
    moments.canonicalMoments().resize(2*nMaxNodes_ - 1);

    // We do not store p0, so canonicalMoments[0] = p1
    scalar p1 = moments.canonicalMoments()[0]; 
    scalar p2 = moments.canonicalMoments()[1];

    scalar alphaCoeff = (1.0 - p1 - 2*p2 + p1*p2)/p2;
    scalar betaCoeff = (p1 - p2 - p1*p2)/p2;

    scalar pJ2n_1 = (betaCoeff + nRegularQuadratureNodes_)
        /(2.0*nRegularQuadratureNodes_ + alphaCoeff + betaCoeff);
    
    scalar pJ2n = nRegularQuadratureNodes_
        /(2.0*nRegularQuadratureNodes_ + 1.0 + alphaCoeff + betaCoeff);

    for
    (
        label i = nRegularQuadratureNodes_; 
        i < nMaxNodes_ && nAdditionalQuadratureNodes_ > 0; 
        i++
    )
    {
        scalar pJ2i_1 = (betaCoeff + i)/(2.0*i + alphaCoeff + betaCoeff);
        scalar pJ2i = i/(2.0*i + 1.0 + alphaCoeff + betaCoeff);

        if (moments.canonicalMoments()[2*nRegularQuadratureNodes_ - 3] <= pJ2i_1 
         || pJ2n_1 >= pJ2i_1)
        {
            moments.canonicalMoments()[2*i - 1] = 
                moments.canonicalMoments()[2*nRegularQuadratureNodes_ - 3]
               *pJ2i_1/pJ2n_1;
        }
        else
        {
            moments.canonicalMoments()[2*i - 1] = 
                (moments.canonicalMoments()[2*nRegularQuadratureNodes_ - 3]
               *(1.0 - pJ2i_1) + pJ2i_1 - pJ2n_1)/(1.0 - pJ2n_1);
        }

        if (moments.canonicalMoments()[2*nRegularQuadratureNodes_ - 2] <= pJ2n
         || pJ2n >= pJ2i)
        {
            moments.canonicalMoments()[2*i] = 
                moments.canonicalMoments()[2*nRegularQuadratureNodes_ - 2]
               *pJ2i/pJ2n;
        }
        else
        {
            moments.canonicalMoments()[2*i] = 
                (moments.canonicalMoments()[2*nRegularQuadratureNodes_ - 2]
               *(1.0 - pJ2i) + pJ2i - pJ2n)/(1.0 - pJ2n);
        }

        moments.zetas()[2*i - 1] = 
            moments.canonicalMoments()[2*i - 1]
           *(1.0 - moments.canonicalMoments()[2*i - 2]);

        moments.zetas()[2*i] = 
            moments.canonicalMoments()[2*i]
           *(1.0 - moments.canonicalMoments()[2*i - 1]);
    }

    alpha[0] = moments.zetas()[0];

    for (label i = 1; i < nMaxNodes_; i++)
    {
        alpha[i] = moments.zetas()[2*i] + moments.zetas()[2*i - 1];
    }

    for (label i = 1; i < nMaxNodes_; i++)
    {
        beta[i] = moments.zetas()[2*i - 1]*moments.zetas()[2*i - 2];
    }
}

// ************************************************************************* //
