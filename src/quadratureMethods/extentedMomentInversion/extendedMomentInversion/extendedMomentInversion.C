/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "extendedMomentInversion.H"
#include "eigenSolver.H"
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
    nMoments_(nMoments),
    nPrimaryNodes_((nMoments_ - 1)/2),
    nSecondaryNodes_(nSecondaryNodes),
    primaryWeights_(nPrimaryNodes_, 0.0),
    primaryAbscissae_(nPrimaryNodes_, 0.0),
    sigma_(0.0),
    secondaryWeights_(nPrimaryNodes_, nSecondaryNodes_),
    secondaryAbscissae_(nPrimaryNodes_, nSecondaryNodes_),
    maxSigmaIter_(dict.lookupOrDefault<label>("maxSigmaIter", 1000)),
    momentsTol_(dict.lookupOrDefault("momentsTol", 1.0e-12)),
    sigmaTol_(dict.lookupOrDefault("sigmaTol", 1.0e-12)),
    targetFunctionTol_(dict.lookupOrDefault("targetFunctionTol", 1.0e-12)),  
    foundUnrealizableSigma_(false),
    nullSigma_(false),
    sigmaBracketed_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::extendedMomentInversion::~extendedMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::extendedMomentInversion::invert(const univariateMomentSet& moments)
{   
    univariateMomentSet m(moments);
    reset();

    // Terminate execution if negative number density is encountered
    if (m[0] < 0.0)
    {
        FatalErrorIn
        (
            "Foam::extendedMomentInversion::invert\n"
            "(\n"
            "   const univariateMomentSet& moments\n"
            ")"
        )   << "The zero-order moment is negative."
            << abort(FatalError);
    }

    // Exclude cases where the zero-order moment is very small to avoid
    // problems in the inversion due to round-off error  
    if (m[0] < SMALL)
    {
        sigma_ = 0.0;
        nullSigma_ = true;

        return;
    }
    
    label nRealizableMoments = m.nRealizableMoments();
   
    if (nRealizableMoments % 2 == 0)
    {
        // If the number of realizable moments is even, we apply the standard
        // QMOM directly to maximize the number of preserved moments.

        // Info << "Even number of realizable moments: using QMOM" << endl;
        // Info << "Moments: " << m << endl;
        // Info << "Invertible: " << m.nInvertibleMoments() << endl;
        // Info << "Realizable: " << m.nRealizableMoments() << endl;
        m.invert();
        sigma_ = 0.0;
        nullSigma_ = true;
        secondaryQuadrature(m);
    }
    else
    {
        // Resizing the moment set to avoid copying again
        m.resize(nRealizableMoments);

        // Local set of starred moments
        univariateMomentSet mStar(nRealizableMoments, 0);

        // Compute target function for sigma = 0
        scalar sigmaLow = 0.0;
        scalar fLow = targetFunction(sigmaLow, m, mStar);
        sigma_ = sigmaLow;
        
        // Check if sigma = 0 is root
        if (mag(fLow) <= targetFunctionTol_)
        {
            m.invert();
            sigma_ = 0.0;
            nullSigma_ = true;
            secondaryQuadrature(m);
            
            return;
        }

        // Compute target function for sigma = sigmaMax 
        scalar sigMax = sigmaMax(m);    
        scalar sigmaHigh = sigMax;
        scalar fHigh = targetFunction(sigmaHigh, m, mStar);

        if (fLow*fHigh > 0)
        {
            // Root not found. Minimize target function in [0, sigma_]
            sigma_ = minimizeTargetFunction(0, sigmaHigh, m, mStar);
            
            // Check if sigma is small and use QMOM
            if (mag(sigma_) < sigmaTol_)
            {
                m.invert();
                sigma_ = 0.0;
                nullSigma_ = true;
                secondaryQuadrature(m);

                return;
            }
            
            targetFunction(sigma_, m, mStar);
            secondaryQuadrature(mStar);

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
                FatalErrorIn
                (
                    "Foam::extendedMomentInversion::invert\n"
                    "(\n"
                    "   const univariateMomentSet& moments\n"
                    ")"
                )<< "Singular value encountered while attempting to find root."
                 << "Moment set = " << m << endl
                 << "sigma = " << sigma_ << endl
                 << "fLow = " << fLow << endl
                 << "fMid = " << fMid << endl
                 << "fHigh = " << fHigh
                 << abort(FatalError);
            }

            sigma_ = sigmaMid + (sigmaMid - sigmaLow)*sign(fLow - fHigh)*fMid/s;
            
            momentsToMomentsStar(sigma_, m, mStar);

            scalar fNew = targetFunction(sigma_, m, mStar);
            scalar dSigma = (sigmaHigh - sigmaLow)/2.0;

            // Check for convergence
            if (mag(fNew) <= targetFunctionTol_ || mag(dSigma) <= sigmaTol_)
            {             
                if (mag(sigma_) < sigmaTol_)
                {
                    m.invert();
                    sigma_ = 0.0;
                    nullSigma_ = true;
                    secondaryQuadrature(m);

                    return;
                }
                
                scalar momentError = normalizedMomentError(sigma_, m, mStar);

                if 
                (
                    momentError < momentsTol_
                )
                {
                    // Found a value of sigma that preserves all the moments
                    secondaryQuadrature(mStar);

                    return;
                }
                else
                {
                    // Root not found. Minimize target function in [0, sigma_]
                    sigma_ = minimizeTargetFunction(0, sigma_, m, mStar);
                    
                    if (mag(sigma_) < sigmaTol_)
                    {
                        m.invert();
                        sigma_ = 0.0;
                        nullSigma_ = true;
                        secondaryQuadrature(m);
                        
                        return;
                    }

                    targetFunction(sigma_, m, mStar);
                    secondaryQuadrature(mStar);

                    return;
                }
            }
            else
            {
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

        FatalErrorIn
        (
            "Foam::extendedMomentInversion::invert\n"
            "(\n"
            "   const univariateMomentSet& moments\n"
            ")"
        )   << "Number of iterations exceeded."
            << abort(FatalError);
    }
}

void Foam::extendedMomentInversion::reset()
{
    foundUnrealizableSigma_ = false;
    nullSigma_ = false;
    sigmaBracketed_ = true;
    
    forAll(primaryWeights_, pNodeI)
    {
        primaryWeights_[pNodeI] = 0.0;
        primaryAbscissae_[pNodeI] = 0.0;
    
        for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
        {
            secondaryWeights_[pNodeI][sNodeI] = 0.0;
            secondaryAbscissae_[pNodeI][sNodeI] = 0.0;
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
        FatalErrorIn
        (
            "Foam::extendedMomentInversion::findExtremumTargetFunction\n"
            "(\n"
            "       const scalar sigmaLow,\n"
            "       const scalar sigmaHigh\n"
            ")"
        )   << "Number of iterations exceeded."
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
    univariateMomentSet approximatedMoments(moments.size(), 0);
    momentsStarToMoments(sigma, approximatedMoments, momentsStar);
    
    //  Info << setprecision (17);
    //  Info << "Approximated moments: " << endl << approximatedMoments;
    //  Info << "Is realizable?" << approximatedMoments.isRealizable() << endl;
    //  Info << "sigma = " << sigma;

    for (label momentI = 0; momentI < moments.size(); momentI++)
    {
        norm += mag(1.0 - approximatedMoments[momentI]/moments[momentI]);
    }

    return sqrt(norm);
}

void Foam::extendedMomentInversion::secondaryQuadrature
(
    const univariateMomentSet& moments
)
{   
    const scalarDiagonalMatrix& pWeights(moments.weights());
    const scalarDiagonalMatrix& pAbscissae(moments.abscissae());
    
    // Copy primary weight and abscissae    
    forAll(pWeights, pNodeI)
    {
        primaryWeights_[pNodeI] = pWeights[pNodeI];
        primaryAbscissae_[pNodeI] = pAbscissae[pNodeI];
    }

    // Set weights and abscissae of unused primary nodes to zero
//     for (label pNodeI = pWeights.size(); pNodeI < nPrimaryNodes_; pNodeI++)
//     {
//         primaryWeights_[pNodeI] = 0.0;
//         primaryAbscissae_[pNodeI] = 0.0;
// 
//         for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
//         {           
//             secondaryWeights_[pNodeI][sNodeI] = 0.0;
//             secondaryAbscissae_[pNodeI][sNodeI] = 0.0;
//         }
//     }
    
    if (!nullSigma_)
    {
        // Coefficients of the recurrence relation
        scalarDiagonalMatrix a(nSecondaryNodes_, 0.0);
        scalarDiagonalMatrix b(nSecondaryNodes_, 0.0);
        
        forAll(pWeights, pNodeI)
        {
            // Compute coefficients of the recurrence relation
            recurrenceRelation(a, b, primaryAbscissae_[pNodeI], sigma_);

            // Define the Jacobi matrix
            scalarSquareMatrix J(nSecondaryNodes_, nSecondaryNodes_, 0.0);

            // Fill diagonal of Jacobi matrix
            forAll(a, aI)
            {
                J[aI][aI] = a[aI];
            }

            // Fill off-diagonal terms of the Jacobi matrix
            for (label bI = 0; bI < nSecondaryNodes_ - 1; bI++)
            {
                J[bI][bI + 1] = Foam::sqrt(b[bI + 1]);
                J[bI + 1][bI] = J[bI][bI + 1];
            }

            // Compute Gaussian quadrature used to find secondary quadrature
            eigenSolver JEig(J, true);

            const scalarDiagonalMatrix& JEigenvaluesRe(JEig.eigenvaluesRe());

            // Compute secondary weights before normalization and calculate sum
            for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
            {
                secondaryWeights_[pNodeI][sNodeI] 
                    = sqr(JEig.eigenvectors()[0][sNodeI]);

                secondaryAbscissae_[pNodeI][sNodeI] = 
                    secondaryAbscissa(primaryAbscissae_[pNodeI], 
                        JEigenvaluesRe[sNodeI], sigma_);
            }
        }
        
        // Set weights and abscissae of unused secondary nodes to zero
        for (label pNodeI = pWeights.size(); pNodeI < nPrimaryNodes_; pNodeI++)
        {
            for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
            {
                secondaryWeights_[pNodeI][sNodeI] = 0.0;
                secondaryAbscissae_[pNodeI][sNodeI] = 0.0;
            }
        }
    }
    else
    {   
        // Manage case with null sigma to avoid redefining source terms
        forAll(pWeights, pNodeI) 
        {
            secondaryWeights_[pNodeI][0] = 1.0;
            secondaryAbscissae_[pNodeI][0] = primaryAbscissae_[pNodeI];

            for (label sNodeI = 1; sNodeI < nSecondaryNodes_; sNodeI++)
            {
                secondaryWeights_[pNodeI][sNodeI] = 0.0;
                secondaryAbscissae_[pNodeI][sNodeI] = 0.0;
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
    momentsStar.invert();
    momentsStar.update();

    scalar lastMoment = moments.last();
    
    return (lastMoment - m2N(sigma, momentsStar))/lastMoment;
}

// ************************************************************************* //
