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
    Copyright (C) 2019-2021 Alberto Passalacqua
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTAbiLITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "extendedMomentInversion.H"
#include "EigenMatrix.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedMomentInversion, 0);
    defineRunTimeSelectionTable(extendedMomentInversion, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedMomentInversion::extendedMomentInversion
(
    const dictionary& dict,
    const label nMoments,
    const label nSecondaryNodes
)
:
    momentInverter_
    (
        univariateMomentInversion::New(dict.subDict("basicQuadrature"))
    ),
    nMoments_(nMoments),
    nPrimaryNodes_((nMoments_ - 1)/2),
    nSecondaryNodes_(nSecondaryNodes),
    primaryWeights_(nPrimaryNodes_, Zero),
    primaryAbscissae_(nPrimaryNodes_, Zero),
    sigma_(0.0),
    secondaryWeights_(nPrimaryNodes_, nSecondaryNodes_),
    secondaryAbscissae_(nPrimaryNodes_, nSecondaryNodes_),
    minMean_(dict.lookupOrDefault<scalar>("minMean", 1.0e-8)),
    minVariance_(dict.lookupOrDefault<scalar>("minVariance", 1.0e-8)),
    maxSigmaIter_(dict.lookupOrDefault<label>("maxSigmaIter", 1000)),
    momentsTol_(dict.lookupOrDefault<scalar>("momentsTol", 1.0e-12)),
    sigmaMin_(dict.lookupOrDefault<scalar>("sigmaMin", 1.0e-6)),
    sigmaTol_(dict.lookupOrDefault<scalar>("sigmaTol", 1.0e-12)),
    targetFunctionTol_
    (
        dict.lookupOrDefault<scalar>("targetFunctionTol", 1.0e-12)
    ),
    foundUnrealizableSigma_(false),
    nullSigma_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedMomentInversion::~extendedMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::extendedMomentInversion::invert(const univariateMomentSet& moments)
{
    univariateMomentSet m(moments);

    reset();

    // Exclude cases where the absolute value of the zero-order moment is very 
    // SMALL to avoid problems in the inversion due to round-off error
    if (mag(m[0]) < SMALL)
    {
        sigma_ = 0.0;
        nullSigma_ = true;

        return;
    }

    // Terminate execution if negative number density is encountered
    if (m[0] < 0.0)
    {
        FatalErrorInFunction
            << "The zero-order moment is negative." << nl
            << "    Moment set: " << m
            << abort(FatalError);
    }

    label nRealizableMoments = m.nRealizableMoments();

    // If the moment set is on the boundary of the moment space, the
    // distribution will be reconstructed by a summation of Dirac delta,
    // and no attempt to use the extended quadrature method of moments is made.
    if (m.isOnMomentSpaceBoundary())
    {
        sigma_ = 0.0;
        nullSigma_ = true;
        momentInverter_().invert(m);

        secondaryQuadrature
        (
            momentInverter_().weights(),
            momentInverter_().abscissae()
        );

        return;
    }

    if (nRealizableMoments % 2 == 0)
    {
        // If the number of realizable moments is even, we apply the standard
        // QMOM directly to maximize the number of preserved moments.
        sigma_ = 0.0;
        nullSigma_ = true;
        momentInverter_().invert(m);

        secondaryQuadrature
        (
            momentInverter_().weights(),
            momentInverter_().abscissae()
        );
    }
    else
    {
        // Do not attempt the EQMOM reconstruction if mean or variance of the
        //  moment set are SMALL to avoid numerical problems. These problems are
        // particularly acute in the calculation of the recurrence relationship
        // of the Jacobi orthogonal polynomials used for the beta kernel density
        // function.
        if (m[1]/m[0] < minMean_ || (m[2]/m[0] - sqr(m[1]/m[0])) < minVariance_)
        {
            sigma_ = 0.0;
            nullSigma_ = true;
            momentInverter_().invert(m);

            secondaryQuadrature
            (
                momentInverter_().weights(),
                momentInverter_().abscissae()
            );

            return;
        }

        // Resizing the moment set to avoid copying again
        m.resize(nRealizableMoments);

        // Local set of starred moments
        univariateMomentSet mStar(nRealizableMoments, m.support());

        // Compute target function for sigma = 0
        scalar sigmaLow = 0.0;
        scalar fLow = targetFunction(sigmaLow, m, mStar);
        sigma_ = sigmaLow;

        // Check if sigma = 0 is root
        if (mag(fLow) <= targetFunctionTol_)
        {
            sigma_ = 0.0;
            nullSigma_ = true;
            momentInverter_().invert(m);

            secondaryQuadrature
            (
                momentInverter_().weights(),
                momentInverter_().abscissae()
            );

            return;
        }

        // Compute target function for sigma = sigmaMax
        scalar sigMax = sigmaMax(m);
        scalar sigmaHigh = sigMax;
        scalar fHigh = targetFunction(sigmaHigh, m, mStar);

        // This should not happen with the new algorithm
        if (fLow*fHigh > 0)
        {
            // Root not found. Minimize target function in [0, sigma_]
            sigma_ = minimizeTargetFunction(0, sigmaHigh, m, mStar);

            // If sigma_ is SMALL, use QMOM
            if (mag(sigma_) < sigmaMin_)
            {
                sigma_ = 0.0;
                nullSigma_ = true;
                momentInverter_().invert(m);

                secondaryQuadrature
                (
                    momentInverter_().weights(),
                    momentInverter_().abscissae()
                );

                return;
            }

            targetFunction(sigma_, m, mStar);
            secondaryQuadrature  // secondary quadrature from mStar
            (
                momentInverter_().weights(),
                momentInverter_().abscissae()
            );

            return;
        }

        // Apply Ridder's algorithm to find sigma
        for (label iter = 0; iter < maxSigmaIter_; iter++)
        {
            scalar fMid, sigmaMid;

            sigmaMid = (sigmaLow + sigmaHigh)/2.0;
            fMid = targetFunction(sigmaMid, m, mStar);

            scalar s = sqrt(sqr(fMid) - fLow*fHigh);

            if (s == 0.0)
            {
                FatalErrorInFunction
                    << "Singular value encountered searching for root." << nl
                    << "    Moment set = " << m << nl
                    << "    sigma = " << sigma_ << nl
                    << "    fLow = " << fLow << nl
                    << "    fMid = " << fMid << nl
                    << "    fHigh = " << fHigh
                    << abort(FatalError);
            }

            sigma_ = sigmaMid + (sigmaMid - sigmaLow)*sign(fLow - fHigh)*fMid/s;

            momentsToMomentsStar(sigma_, m, mStar);

            scalar fNew = targetFunction(sigma_, m, mStar);
            scalar dSigma = (sigmaHigh - sigmaLow)/2.0;

            // Check for convergence
            if (mag(fNew) <= targetFunctionTol_ || mag(dSigma) <= sigmaTol_)
            {
                // Root finding converged

                // If sigma_ is SMALL, use QMOM
                if (mag(sigma_) < sigmaMin_)
                {
                    sigma_ = 0.0;
                    nullSigma_ = true;
                    momentInverter_().invert(m);

                    secondaryQuadrature
                    (
                        momentInverter_().weights(),
                        momentInverter_().abscissae()
                    );

                    return;
                }

                scalar momentError = normalizedMomentError(sigma_, m, mStar);

                if
                (
                    momentError < momentsTol_
                )
                {
                    // Found a value of sigma that preserves all the moments
                    secondaryQuadrature  // Secondary quadrature from mStar
                    (
                        momentInverter_().weights(),
                        momentInverter_().abscissae()
                    );

                    return;
                }
                else
                {
                    // Root not found. Minimize target function in [0, sigma_]
                    sigma_ = minimizeTargetFunction(0, sigma_, m, mStar);

                    // If sigma_ is SMALL, use QMOM
                    if (mag(sigma_) < sigmaMin_)
                    {
                        sigma_ = 0.0;
                        nullSigma_ = true;
                        momentInverter_().invert(m);

                        secondaryQuadrature
                        (
                            momentInverter_().weights(),
                            momentInverter_().abscissae()
                        );

                        return;
                    }

                    targetFunction(sigma_, m, mStar);

                    secondaryQuadrature // Secondary quadrature from  mStar
                    (
                        momentInverter_().weights(),
                        momentInverter_().abscissae()
                    );

                    return;
                }
            }
            else
            {
                // Root finding did not converge. Refine search.

                if (fNew*fMid < 0 && sigma_ < sigmaMid)
                {
                    sigmaLow = sigma_;
                    fLow = fNew;
                    sigmaHigh = sigmaMid;
                    fHigh = fMid;
                }
                else if (fNew*fMid < 0 && sigma_ > sigmaMid)
                {
                    sigmaLow = sigmaMid;
                    fLow = fMid;
                    sigmaHigh = sigma_;
                    fHigh = fNew;
                }
                else if (fNew*fLow < 0)
                {
                    sigmaHigh = sigma_;
                    fHigh = fNew;
                }
                else if (fNew*fHigh < 0)
                {
                    sigmaLow = sigma_;
                    fLow = fNew;
                }
            }
        }

        FatalErrorInFunction
            << "Number of iterations exceeded." << nl
            << "    Max allowed iterations = " << maxSigmaIter_
            << abort(FatalError);
    }
}

void Foam::extendedMomentInversion::reset()
{
    foundUnrealizableSigma_ = false;
    nullSigma_ = false;

    forAll(primaryWeights_, pNodei)
    {
        primaryWeights_[pNodei] = 0.0;
        primaryAbscissae_[pNodei] = 0.0;

        for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
        {
            secondaryWeights_[pNodei][sNodei] = 0.0;
            secondaryAbscissae_[pNodei][sNodei] = 0.0;
        }
    }
}

Foam::scalar Foam::extendedMomentInversion::minimizeTargetFunction
(
    scalar sigmaLow,
    scalar sigmaHigh,
    const univariateMomentSet& moments,
    univariateMomentSet& momentsStar
)
{
    const scalar goldenRatio = (sqrt(5.0) - 1.0)/2.0;

    scalar a = sigmaLow;
    scalar b = sigmaHigh;
    scalar x = b - goldenRatio*(b - a);
    scalar y = a + goldenRatio*(b - a);

    label iter = 0;

    while (mag (x - y) > sigmaTol_ && iter < maxSigmaIter_)
    {
        // Square the target function to find closest value to zero,
        // independently from the sign of the function
        scalar fx = sqr(targetFunction(x, moments, momentsStar));
        scalar fy = sqr(targetFunction(y, moments, momentsStar));

        if (fx < fy)
        {
            b = y;
            y = x;
            x = b - goldenRatio*(b - a);
        }
        else
        {
            a = x;
            x = y;
            y = a + goldenRatio*(b - a);
        }

        iter++;
    }

    if (iter > maxSigmaIter_)
    {
        FatalErrorInFunction
            << "Number of iterations exceeded." << nl
            << "    Max allowed iterations = " << maxSigmaIter_
            << abort(FatalError);
    }

    return (a + b)/2.0;
}

Foam::scalar Foam::extendedMomentInversion::normalizedMomentError
(
    scalar sigma,
    const univariateMomentSet& moments,
    univariateMomentSet& momentsStar
)
{
    scalar norm = 0.0;

    targetFunction(sigma, moments, momentsStar);

    univariateMomentSet approximatedMoments(moments.size(), moments.support());

    momentsStarToMoments(sigma, approximatedMoments, momentsStar);

    for (label momenti = 0; momenti < moments.size(); momenti++)
    {
        norm += mag(1.0 - approximatedMoments[momenti]/moments[momenti]);
    }

    return sqrt(norm);
}

void Foam::extendedMomentInversion::secondaryQuadrature
(
    const scalarList& pWeights,
    const scalarList& pAbscissae
)
{
    // Copy primary weights and abscissae
    forAll(pWeights, pNodei)
    {
        primaryWeights_[pNodei] = pWeights[pNodei];
        primaryAbscissae_[pNodei] = pAbscissae[pNodei];
    }

    if (!nullSigma_)
    {
        // Coefficients of the recurrence relation
        scalarDiagonalMatrix a(nSecondaryNodes_, Zero);
        scalarDiagonalMatrix b(nSecondaryNodes_, Zero);

        forAll(pWeights, pNodei)
        {
            // Compute coefficients of the recurrence relation
            recurrenceRelation(a, b, primaryAbscissae_[pNodei], sigma_);

            // Define the Jacobi matrix
            scalarSquareMatrix J(nSecondaryNodes_, Zero);

            // Fill diagonal of Jacobi matrix
            forAll(a, ai)
            {
                J[ai][ai] = a[ai];
            }

            // Fill off-diagonal terms of the Jacobi matrix
            for (label bi = 0; bi < nSecondaryNodes_ - 1; bi++)
            {
                J[bi][bi + 1] = Foam::sqrt(b[bi + 1]);
                J[bi + 1][bi] = J[bi][bi + 1];
            }

            // Compute Gaussian quadrature used to find secondary quadrature
            EigenMatrix<scalar> JEig(J, true);

            const scalarDiagonalMatrix& JEigenvaluesRe(JEig.EValsRe());

            // Compute secondary weights before normalization and calculate sum
            for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
            {
                secondaryWeights_[pNodei][sNodei] 
                    = sqr(JEig.EVecs()[0][sNodei]);

                secondaryAbscissae_[pNodei][sNodei] =
                    secondaryAbscissa(primaryAbscissae_[pNodei],
                        JEigenvaluesRe[sNodei], sigma_);
            }
        }

        // Set weights and abscissae of unused nodes to zero
        for (label pNodei = pWeights.size(); pNodei < nPrimaryNodes_; pNodei++)
        {
            primaryWeights_[pNodei] = 0.0;
            primaryAbscissae_[pNodei] = 0.0;

            for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
            {
                secondaryWeights_[pNodei][sNodei] = 0.0;
                secondaryAbscissae_[pNodei][sNodei] = 0.0;
            }
        }
    }
    else
    {
        // Manage case with null sigma to avoid redefining source terms
        forAll(pWeights, pNodei)
        {
            secondaryWeights_[pNodei][0] = 1.0;
            secondaryAbscissae_[pNodei][0] = primaryAbscissae_[pNodei];

            for (label sNodei = 1; sNodei < nSecondaryNodes_; sNodei++)
            {
                secondaryWeights_[pNodei][sNodei] = 0.0;
                secondaryAbscissae_[pNodei][sNodei] = 0.0;
            }
        }

        // Set weights and abscissae of unused nodes to zero
        for (label pNodei = pWeights.size(); pNodei < nPrimaryNodes_; pNodei++)
        {
            primaryWeights_[pNodei] = 0.0;
            primaryAbscissae_[pNodei] = 0.0;

            for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
            {
                secondaryWeights_[pNodei][sNodei] = 0.0;
                secondaryAbscissae_[pNodei][sNodei] = 0.0;
            }
        }
    }
}

Foam::scalar Foam::extendedMomentInversion::targetFunction
(
    scalar sigma,
    const univariateMomentSet& moments,
    univariateMomentSet& momentsStar
)
{
    momentsToMomentsStar(sigma, moments, momentsStar);

    momentInverter_().invert(momentsStar);
    momentsStar.update
    (
        momentInverter_().weights(),
        momentInverter_().abscissae()
    );

    scalar lastMoment = moments.last();

    return (lastMoment - m2N(sigma, momentsStar))/lastMoment;
}

// ************************************************************************* //
