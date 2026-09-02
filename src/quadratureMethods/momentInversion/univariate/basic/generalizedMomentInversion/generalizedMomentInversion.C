/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2026 Alberto Passalacqua
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
#include "supportType.H"

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
            << "The type of NDF for RPlus must be gamma or" << nl
            << "lognormal. The current value is " << ndfTypeRPlus_
            << abort(FatalError);
    }
}


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
    const supportType& support = moments.support();

    #ifdef FULLDEBUG
        Info << "Support = " << supportTypeToWord(support) << endl;
    #endif

    if (support == supportType::R)
    {
        correctRecurrenceR(alpha, beta);
    }
    else if (support == supportType::RPlus)
    {
        correctRecurrenceRPlus(moments, alpha, beta);
    }
    else if (support == supportType::ZeroOne)
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

    // Resize list of weights and abscissae.
    // Note: the lists for the alpha and beta coefficients of the recurrence
    //       relationship, and those for the zeta_k and the canonical moments,
    //       do NOT need to be resized: univariateMomentSet allocates them
    //       with room for the additional quadrature nodes.
    weights_.setSize(nMaxNodes_);
    abscissae_.setSize(nMaxNodes_);

    if (moments.zetas().size() < 2*nMaxNodes_ - 1)
    {
        FatalErrorInFunction
            << "The moment set does not reserve room for the zeta_k of the "
            << "additional quadrature nodes." << nl
            << "    Number of quadrature nodes: " << nMaxNodes_ << nl
            << "    Size of the zeta list: " << moments.zetas().size() << nl
            << "    Required size: " << 2*nMaxNodes_ - 1 << nl
            << nl
            << "The moment set was built for a quadrature with fewer nodes "
            << "than GQMOM requires. Check that the number of additional "
            << "quadrature points of the fieldMomentInversion that owns this "
            << "moment set was derived from GQMOM, and not from the "
            << "quadrature selected in another dictionary of the case." << nl
            << exit(FatalError);
    }

    #ifdef FULLDEBUG
        Info << "nMaxNodes = " << nMaxNodes_ << endl
            << "nRegularQuadratureNodes = "
            << nRegularQuadratureNodes_ << endl
            << "nAdditionalQuadratureNodes = "
            << nAdditionalQuadratureNodes_ << endl;
    #endif
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

    #ifdef FULLDEBUG
        Info << "Corrected alpha: " << alpha << endl;
        Info << "Corrected beta: " << beta << endl;
    #endif
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

    // Take a reference to zetas and use it instead than
    // accessing moments.zetas() directly.
    scalarList& zetas(moments.zetas());

    if (ndfTypeRPlus_ == "gamma")
    {
        // The moments are read through the const accessor, so that the
        // realizability check, and with it the zeta chain extended below, is
        // not invalidated
        const scalar m1sqr = sqr(moments.moment(1));

        const scalar alphaCoeff =
            m1sqr/(moments.moment(0)*moments.moment(2) - m1sqr) - 1.0;

        for
        (
            label i = nRegularQuadratureNodes_;
            i < nMaxNodes_ && nAdditionalQuadratureNodes_ > 0;
            i++
        )
        {
            zetas[2*i - 1] =
                (i + alphaCoeff)*zetas[2*nRegularQuadratureNodes_ - 3]
               /(nRegularQuadratureNodes_ - 1 + alphaCoeff);

            zetas[2*i] =
                (i + 1)*zetas[2*nRegularQuadratureNodes_ - 2]
               /(nRegularQuadratureNodes_);

            #ifdef FULLDEBUG
                Info << "zetas[2*i-1] = " << zetas[2*i - 1] << endl;
                Info << "2i-1 = " << 2*i - 1 << endl;
                Info << "zetas[2*i] = " << zetas[2*i] << endl;
                Info << "2i = " << 2*i << endl;
            #endif
        }
    }
    else if (ndfTypeRPlus_ == "lognormal")
    {
        const scalar eta =
            sqrt
            (
                moments.moment(0)*moments.moment(2)/sqr(moments.moment(1))
            );

        for
        (
            label i = nRegularQuadratureNodes_;
            i < nMaxNodes_ && nAdditionalQuadratureNodes_ > 0;
            i++
        )
        {
            zetas[2*i - 1] =
                pow(eta, 2*(i + 1 - nRegularQuadratureNodes_))
               *(
                    (pow(eta, 2*(i+1)) - 1.0)
                   /(pow(eta, 2*nRegularQuadratureNodes_) - 1.0)
                )
               *zetas[2*nRegularQuadratureNodes_ - 3];

            zetas[2*i] =
                pow(eta, 4*(i + 1 - nRegularQuadratureNodes_))
               *zetas[2*nRegularQuadratureNodes_ - 2];

            #ifdef FULLDEBUG
                Info << "zetas[2*i-1] = " << zetas[2*i - 1] << endl;
                Info << "2i-1 = " << 2*i - 1 << endl;
                Info << "zetas[2*i] = " << zetas[2*i] << endl;
                Info << "2i = " << 2*i << endl;
            #endif
        }
    }

    // The zeta chain is the primary representation of the recurrence
    // relationship for measures with support over R+ and [0, 1]: recover
    // alpha and beta from the chain extended above
    moments.zetasToRecurrence(2*nMaxNodes_ - 1, alpha, beta);

    #ifdef FULLDEBUG
        Info << "Corrected alpha: " << alpha << endl;
        Info << "Corrected beta: " << beta << endl;
    #endif
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

    scalarList& zetas(moments.zetas());
    scalarList& canonicalMoments(moments.canonicalMoments());

    // We do not store p0, so canonicalMoments[0] = p1
    scalar p1 = canonicalMoments[0];
    scalar p2 = canonicalMoments[1];

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

        if (canonicalMoments[2*nRegularQuadratureNodes_ - 3] <= pJ2n_1
         || pJ2n_1 >= pJ2i_1)
        {
            canonicalMoments[2*i - 1] =
                canonicalMoments[2*nRegularQuadratureNodes_ - 3]
               *pJ2i_1/pJ2n_1;
        }
        else
        {
            canonicalMoments[2*i - 1] =
                (canonicalMoments[2*nRegularQuadratureNodes_ - 3]
               *(1.0 - pJ2i_1) + pJ2i_1 - pJ2n_1)/(1.0 - pJ2n_1);
        }

        if (canonicalMoments[2*nRegularQuadratureNodes_ - 2] <= pJ2n
         || pJ2n >= pJ2i)
        {
            canonicalMoments[2*i] =
                canonicalMoments[2*nRegularQuadratureNodes_ - 2]
               *pJ2i/pJ2n;
        }
        else
        {
            canonicalMoments[2*i] =
                (canonicalMoments[2*nRegularQuadratureNodes_ - 2]
               *(1.0 - pJ2i) + pJ2i - pJ2n)/(1.0 - pJ2n);
        }

        zetas[2*i - 1] =
            canonicalMoments[2*i - 1]
           *(1.0 - canonicalMoments[2*i - 2]);

        zetas[2*i] =
            canonicalMoments[2*i]
           *(1.0 - canonicalMoments[2*i - 1]);

        #ifdef FULLDEBUG
            Info << "canonicalMoments[2*i-1] = " << canonicalMoments[2*i - 1] << endl;
            Info << "zetas[2*i-1] = " << zetas[2*i - 1] << endl;
            Info << "2i-1 = " << 2*i - 1 << endl;
            Info << "canonicalMoments[2*i] = " << canonicalMoments[2*i] << endl;
            Info << "zetas[2*i] = " << zetas[2*i] << endl;
            Info << "2i = " << 2*i << endl;
        #endif
    }

    // The zeta chain is the primary representation of the recurrence
    // relationship for measures with support over R+ and [0, 1]: recover
    // alpha and beta from the chain extended above
    moments.zetasToRecurrence(2*nMaxNodes_ - 1, alpha, beta);

    #ifdef FULLDEBUG
        Info << "Corrected alpha: " << alpha << endl;
        Info << "Corrected beta: " << beta << endl;
    #endif
}

// ************************************************************************* //
