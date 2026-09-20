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
    nMaxNodes_(nMaxNodes),
    nRealizableMoments_(0)
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
    scalarList& beta
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

bool Foam::generalizedMomentInversion::canExtendRecurrence
(
    const univariateMomentSet& moments
) const
{
    // The model number density functions that extend the recurrence
    // relationship to the additional nodes are singular for a moment vector
    // without spread: the gamma coefficient divides by m0 m2 - m1^2, the
    // lognormal one by eta^(2n) - 1 with eta^2 = m0 m2/m1^2, which underflows
    // to zero for the same reason, and the Jacobi one by the second canonical
    // moment. A moment vector of that kind is a Dirac delta, which the
    // regular quadrature nodes already represent exactly, so there is nothing
    // for the additional nodes to add and Gauss is used instead.
    const supportType& support = moments.support();

    if (support == supportType::RPlus)
    {
        const scalar m1 = moments.moment(1);
        const scalar variance = moments.moment(0)*moments.moment(2) - sqr(m1);

        return m1 > 0 && variance > SMALL*sqr(m1);
    }

    if (support == supportType::ZeroOne)
    {
        return moments.canonicalMoments()[1] > SMALL;
    }

    // Over R the recurrence relationship is extended from the mean of the
    // alpha coefficients and a power law on the beta ones, which is regular
    // for every moment vector with more than three realizable moments
    return true;
}


Foam::label Foam::generalizedMomentInversion::nAdditionalQuadraturePoints
(
    const label nMoments
) const
{
    // The moments determine at most nMoments/2 quadrature nodes, for an even
    // and for an odd number of moments alike, and the recurrence relationship
    // is extended with a model number density function to the remaining ones
    const label nAdditionalNodes = nMaxNodes_ - nMoments/2;

    if (nAdditionalNodes < 0)
    {
        FatalErrorInFunction
            << "Generalized QMOM is asked for fewer quadrature nodes than "
            << "the moments determine." << nl
            << "    Number of moments: " << nMoments << nl
            << "    Number of quadrature nodes: " << nMaxNodes_ << nl
            << "    Number of nodes of regular QMOM: " << nMoments/2 << nl
            << nl
            << "The moments beyond the ones the requested nodes support "
            << "would be lost, and the quadrature would not conserve them. "
            << "The number of nodes cannot be raised to the one of regular "
            << "QMOM either, because the nodes in excess of the ones the "
            << "case declares are discarded when the quadrature is stored."
            << nl
            << "Declare at least " << nMoments/2 << " nodes, or reduce the "
            << "number of moments." << nl
            << exit(FatalError);
    }

    return nAdditionalNodes;
}


void Foam::generalizedMomentInversion::calcNQuadratureNodes
(
    univariateMomentSet& moments
)
{
    // Trigger calculations of zeta_k by computing the number of realizable
    // moments.
    label nRealizableMoments = moments.nRealizableMoments();

    nRealizableMoments_ = nRealizableMoments;

    nRegularQuadratureNodes_
        = (nRealizableMoments % 2 != 0)
        ? label((nRealizableMoments - 1)/2.0)
        : label(nRealizableMoments/2.0);

    // Note: nRegularQuadratureNodes_ is at most nMoments/2, because the
    //       realizable moments are at most all of them, and the owner of the
    //       moment set has checked that nMaxNodes_ is not smaller than that,
    //       so the number of additional nodes below cannot be negative.
    if (nRealizableMoments > 3 && canExtendRecurrence(moments))
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

    if
    (
        nAdditionalQuadratureNodes_ > 0
     && moments.zetas().size() < 2*nMaxNodes_ - 1
    )
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
    univariateMomentSet& moments
)
{
    (*this).univariateMomentInversion::invert
        (
            moments
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

    // The a_i the moments determine are as many as the regular nodes, and
    // the ones beyond them are their mean
    scalar an = 0;

    for (label i = 0; i < nRegularQuadratureNodes_; i++)
    {
        an += alpha[i];
    }

    an /= nRegularQuadratureNodes_;

    for (label i = nRegularQuadratureNodes_; i < nNodes_; i++)
    {
        alpha[i] = an;
    }

    // The last b_i the moments determine: b_(n-1) for 2n moments and b_n
    // for 2n + 1, which the moment of order 2n gives. The ones beyond it
    // follow the power law of the Hermite polynomials from it, so that the
    // moment of order 2n is kept where it is known.
    const label lastB = (nRealizableMoments_ - 1)/2;

    for (label i = lastB + 1; i < nNodes_; i++)
    {
        beta[i] = beta[lastB]*pow(scalar(i)/scalar(lastB), nu_);
    }

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

    // The zeta_k are stored from zeta_1, so zeta_j is in zetas[j - 1]. The
    // realizable moments determine zeta_1 ... zeta_nKnown. The chain of the
    // odd orders, zeta_(2k-1), and the one of the even orders, zeta_(2k),
    // each continue from the last value the moments determine, with the
    // ratio the model distribution gives between the index k of the value
    // sought and the one of that value.
    const label nKnown = nRealizableMoments_ - 1;

    // Index k of the last known zeta_(2k-1) and zeta_(2k)
    const label kOdd = (nKnown + 1)/2;
    const label kEven = nKnown/2;

    const scalar lastOdd = zetas[2*kOdd - 2];
    const scalar lastEven = zetas[2*kEven - 1];

    if (ndfTypeRPlus_ == "gamma")
    {
        // The moments are read through the const accessor, so that the
        // realizability check, and with it the zeta chain extended below, is
        // not invalidated
        const scalar m1sqr = sqr(moments.moment(1));

        const scalar alphaCoeff =
            m1sqr/(moments.moment(0)*moments.moment(2) - m1sqr) - 1.0;

        // Generalized Laguerre polynomials:
        //     zeta_(2k-1) = (k + alpha)/beta, zeta_(2k) = k/beta
        for (label j = nKnown + 1; j <= 2*nMaxNodes_ - 1; j++)
        {
            const label k = (j + 1)/2;

            zetas[j - 1] =
                (j % 2 == 1)
              ? lastOdd*(k + alphaCoeff)/(kOdd + alphaCoeff)
              : lastEven*scalar(k)/scalar(kEven);
        }
    }
    else if (ndfTypeRPlus_ == "lognormal")
    {
        const scalar eta =
            sqrt
            (
                moments.moment(0)*moments.moment(2)/sqr(moments.moment(1))
            );

        // Stieltjes-Wigert polynomials:
        //     zeta_(2k-1) = exp(mu) eta^(4k-3),
        //     zeta_(2k) = exp(mu) eta^(2k-1) (eta^(2k) - 1)
        for (label j = nKnown + 1; j <= 2*nMaxNodes_ - 1; j++)
        {
            const label k = (j + 1)/2;

            zetas[j - 1] =
                (j % 2 == 1)
              ? lastOdd*pow(eta, scalar(4*(k - kOdd)))
              : lastEven*pow(eta, scalar(2*(k - kEven)))
               *(pow(eta, scalar(2*k)) - 1.0)
               /(pow(eta, scalar(2*kEven)) - 1.0);
        }
    }

    #ifdef FULLDEBUG
        Info << "Extended zetas: " << zetas << endl;
    #endif

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

    // We do not store p0, so p_j is in canonicalMoments[j - 1], as zeta_j
    // is in zetas[j - 1]
    const scalar p1 = canonicalMoments[0];
    const scalar p2 = canonicalMoments[1];

    const scalar alphaCoeff = (1.0 - p1 - 2*p2 + p1*p2)/p2;
    const scalar betaCoeff = (p1 - p2 - p1*p2)/p2;

    // Canonical moments of the Jacobi polynomials:
    //     p_(2k-1) = (beta + k)/(2k + alpha + beta),
    //     p_(2k) = k/(2k + 1 + alpha + beta)
    auto pJ = [&](const label j)
    {
        const label k = (j + 1)/2;

        return
            (j % 2 == 1)
          ? (betaCoeff + k)/(2.0*k + alphaCoeff + betaCoeff)
          : scalar(k)/(2.0*k + 1.0 + alphaCoeff + betaCoeff);
    };

    // The realizable moments determine p_1 ... p_nKnown. The chain of the
    // odd orders and the one of the even orders each continue from the last
    // value the moments determine.
    const label nKnown = nRealizableMoments_ - 1;

    const label jOdd = (nKnown % 2 == 1) ? nKnown : nKnown - 1;
    const label jEven = (nKnown % 2 == 0) ? nKnown : nKnown - 1;

    for (label j = nKnown + 1; j <= 2*nMaxNodes_ - 1; j++)
    {
        const label jLast = (j % 2 == 1) ? jOdd : jEven;

        const scalar pLast = canonicalMoments[jLast - 1];
        const scalar pJLast = pJ(jLast);
        const scalar pJj = pJ(j);

        // The value sought scales with the one of the Jacobi polynomials
        // from zero where the last known one is below its own, and from one
        // where it is above, so that it stays inside (0, 1)
        canonicalMoments[j - 1] =
            (pLast <= pJLast || pJLast >= pJj)
          ? pLast*pJj/pJLast
          : (pLast*(1.0 - pJj) + pJj - pJLast)/(1.0 - pJLast);

        zetas[j - 1] =
            canonicalMoments[j - 1]*(1.0 - canonicalMoments[j - 2]);
    }

    #ifdef FULLDEBUG
        Info << "Extended canonical moments: " << canonicalMoments << endl;
        Info << "Extended zetas: " << zetas << endl;
    #endif

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
