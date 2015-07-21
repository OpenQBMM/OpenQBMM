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
    lastMomentI_(nMoments_ - 1),
    nPrimaryNodes_((nMoments_ - 1)/2),
    nSecondaryNodes_(nSecondaryNodes),
    approximatedMoments_(nMoments_, 0.0),
    momentsStar_(nMoments_, 0.0),
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
    
    foundUnrealizableSigma_ = false;
    nullSigma_ = false;
    sigmaBracketed_ = true;
    
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
        return;
    }
    
    label nRealizableMoments = m.nRealizableMoments();
    
    if (nRealizableMoments % 2 == 0)
    {
        // If the number of realizable moments is even, we apply the standard
        // QMOM directly to maximize the number of preserved moments.

        Info << "Even number of realizable moments: using QMOM" << endl;
        Info << "Moments: " << m << endl;
        Info << "Invertible: " << m.nInvertibleMoments() << endl;
        Info << "Realizable: " << m.nRealizableMoments() << endl;
        m.invert();
        nullSigma_ = true;
        secondaryQuadrature(m);
    }
    else
    {
        scalar sigmaLow = 0.0;
        scalar fLow = targetFunction(sigmaLow, m);
        sigma_ = sigmaLow;
        
//         if (fLow < targetFunctionTol_)
//         {
//             nullSigma_ = true;
//             // TODO Transfer w/a        
//             return;
//         }
        
        // Find maximum value of sigma and check if it is root
        scalar sigMax = sigmaMax(m);    
        scalar sigmaHigh = sigMax;
        scalar fHigh = targetFunction(sigmaHigh, m);
        
        scalar momentErrorSigmaMax = normalizedMomentError(sigMax, m);
        
        if (fLow*fHigh > 0)
        {
            // Attempting to bracket
            scalar sigmaExtremumTargetFunction 
                    = findExtremumTargetFunction(sigmaLow, sigmaHigh, m);

            scalar fMin = targetFunction(sigmaExtremumTargetFunction, m);

            if (fMin*fLow < 0)
            {
                fLow = targetFunction(sigmaLow, m);
                sigmaHigh = sigmaExtremumTargetFunction;
                fHigh = targetFunction(sigmaHigh, m);
            }
            else
            {
                sigmaBracketed_ = false;

                sigma_ = sigmaExtremumTargetFunction;

                momentsToMomentsStar(sigma_, m);
                momentsStar_.invert();
                scalar momentError = normalizedMomentError(sigma_, m);
                
                // Use sigma from minimization
                
                if (momentError > momentsTol_)
                {
                    Info<< "The moment set is not preserved.\n" 
                        << "sigma = " << sigma_ << endl
                        << "Moment error = " << momentError << endl;
                }
                
                secondaryQuadrature(momentsStar_);

                return;
            }
        }
        
        for (label iter = 0; iter < maxSigmaIter_; iter++)
        {
            scalar sigmaMid = (sigmaLow + sigmaHigh)/2.0;
            scalar fMid = targetFunction(sigmaMid, m);

    //         Info << "sigmaMid = " << sigmaMid << endl;
    //         Info << "sigmaLow = " << sigmaLow << " fLow = " << fLow << endl;
    //         Info << "fMid = " << fMid << endl;
    //         Info << "sigmaHigh = " << sigmaHigh << "fHigh = " << fHigh << endl;
        
            scalar s = sqrt(sqr(fMid) - fLow*fHigh);

            if (s == 0.0)
            {
                FatalErrorIn
                (
                    "Foam::extendedMomentInversion::invert\n"
                    "(\n"
                    "   const univariateMomentSet& moments\n"
                    ")"
                )   << "Singular value encountered while attempting to find root."
                    << abort(FatalError);
            }

            sigma_ = sigmaMid + (sigmaMid - sigmaLow)*sign(fLow - fHigh)*fMid/s;
            
            // NOTE Make momentStar_ local?
            momentsToMomentsStar(sigma_, m);

            if (!momentsStar_.isRealizable() && !foundUnrealizableSigma_)
            {
                foundUnrealizableSigma_ = true;
            }

            scalar fNew = targetFunction(sigma_, m);
            scalar dSigma = (sigmaHigh - sigmaLow)/2.0;

            // Check for convergence
            if (mag(fNew) <= targetFunctionTol_ || mag(dSigma) <= sigmaTol_)
            {
                scalar momentError = normalizedMomentError(sigma_, m);

                if 
                (
                    foundUnrealizableSigma_
                    && !momentsStar_.isRealizable()
                    && momentError > momentErrorSigmaMax
                )
                {
                    // If a value of sigma that leads to unrealizable moments is
                    // found, use the value of sigma that minimizes the error 
                    // on the moments TODO Check
                    sigma_ = sigMax;
                    momentsToMomentsStar(sigma_, m);
                    momentsStar_.invert();
                    momentError = normalizedMomentError(sigma_, m);
                }

                if (momentError > momentsTol_)
                {
                    Info<< "The moment set is not preserved." << endl;
                    Info<< "Moment error = " << momentError << endl;
                    Info<< "Sigma = " << sigma_ << endl;
                    Info<< "Target function = " << fNew << endl;
                    Info<< "Sigma change = " << dSigma << endl;
                    
                    // Checking if there is another root
                    sigmaLow = 0.0;
                    fLow = targetFunction(sigmaLow, m);
                    
                    scalar sigmaExtremumTargetFunction 
                        = findExtremumTargetFunction(sigmaLow, sigma_, m);
                        
                    scalar fMin = targetFunction(sigmaExtremumTargetFunction, m);
                    
                    Info << "sigmaMinNew = " << sigmaExtremumTargetFunction 
                         << endl << "fMinNew = " << fMin << endl;
                    
                    if (fMin*fLow <= 0.0)
                    {
                        Info << "Finding other root" << endl;

                        fHigh = fMin;
                        sigmaHigh = sigmaExtremumTargetFunction;
                        iter = 0;
                    }
                    else
                    {
                        Info << "Here!" << endl;
                        secondaryQuadrature(momentsStar_);
                        return;
                    }
                }
                else
                {
                    secondaryQuadrature(momentsStar_);
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
    
//     forAll(primaryWeights_, wI)
//     {
//         primaryWeights_[wI] = 0.0;
//         primaryAbscissae_[wI] = 0.0;
//     }
}

Foam::scalar Foam::extendedMomentInversion::findExtremumTargetFunction
(
    scalar sigmaLow,
    scalar sigmaHigh,
    const univariateMomentSet& moments
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
        scalar fx = sqr(targetFunction(x, moments));
        scalar fy = sqr(targetFunction(y, moments));
    
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
    const univariateMomentSet& moments
)
{
    scalar norm = 0.0;

    targetFunction(sigma, moments);
    
    //Info << setprecision (17);
    //Info << "Approximated moments: " << endl << approximatedMoments_;
    //Info << "Is realizable?" << approximatedMoments_.isRealizable() << endl;

    for (label momentI = 0; momentI < nMoments_; momentI++)
    {
        norm += magSqr((moments[momentI] - approximatedMoments_[momentI])
                /moments[momentI]);
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
        
    forAll(pWeights, pNodeI)
    {
        primaryWeights_[pNodeI] = pWeights[pNodeI];
        primaryAbscissae_[pNodeI] = pAbscissae[pNodeI];
    }
    
    if (!nullSigma_)
    {
        // Coefficients of the recurrence relation
        scalarDiagonalMatrix a(nSecondaryNodes_, 0.0);
        scalarDiagonalMatrix b(nSecondaryNodes_, 0.0);
        
        forAll(primaryWeights_, pNodeI)
        {
            // Compute coefficients of the recurrence relation
            recurrenceRelation(a, b, primaryAbscissae_[pNodeI]);

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
    }
    else
    {       
        // Manage case with null sigma to avoid redefining source terms
        forAll(primaryWeights_, pNodeI)
        {           
            if (primaryWeights_[pNodeI] > 0)
            {
                secondaryWeights_[pNodeI][0] = 1.0;
                secondaryAbscissae_[pNodeI][0] = primaryAbscissae_[pNodeI];
            }
            else
            {
                secondaryWeights_[pNodeI][0] = 0.0;
                secondaryAbscissae_[pNodeI][0] = 0.0;
            }

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
    const univariateMomentSet& moments
)
{
    momentsToMomentsStar(sigma, moments);
    momentsStar_.invert();
    momentsStar_.update();

    // NOTE: make approximatedMoments local?
    momentsStarToMoments(sigma, approximatedMoments_);

    return (moments[lastMomentI_] - approximatedMoments_[lastMomentI_])
          /moments[lastMomentI_];
}

// ************************************************************************* //
