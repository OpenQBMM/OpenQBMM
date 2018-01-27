/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

#include "hyperbolicConditionalMomentInversion.H"

#include "mappedLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::labelListList
Foam::hyperbolicConditionalMomentInversion::hyperbolicConditionalMomentInversion
::twoDimMomentOrders_ =
{
    {0, 0},
    {1, 0},
    {0, 1},
    {2, 0},
    {1, 1},
    {0, 2},
    {3, 0},
    {2, 1},
    {1, 2},
    {0, 3},
    {4, 0},
    {0, 4}
};

const Foam::labelListList
Foam::hyperbolicConditionalMomentInversion::hyperbolicConditionalMomentInversion
::threeDimMomentOrders_ =
{
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1},
    {2, 0, 0},
    {1, 1, 0},
    {1, 0, 1},
    {0, 2, 0},
    {0, 1, 1},
    {0, 0, 2},
    {3, 0, 0},
    {2, 1, 0},
    {2, 0, 1},
    {1, 2, 0},
    {1, 1, 1},
    {1, 0, 2},
    {0, 3, 0},
    {0, 2, 1},
    {0, 1, 2},
    {0, 0, 3},
    {4, 0, 0},
    {0, 4, 0},
    {0, 0, 4}
};


const Foam::labelListList
Foam::hyperbolicConditionalMomentInversion::hyperbolicConditionalMomentInversion
::twoDimNodeIndexes_ =
{
    {1, 1},
    {1, 2},
    {1, 3},
    {2, 1},
    {2, 2},
    {2, 3},
    {3, 1},
    {3, 2},
    {3, 3}
};

const Foam::labelListList
Foam::hyperbolicConditionalMomentInversion::hyperbolicConditionalMomentInversion
::threeDimNodeIndexes_ =
{
    {1, 1, 1},
    {1, 1, 2},
    {1, 1, 3},
    {1, 2, 1},
    {1, 2, 2},
    {1, 2, 3},
    {1, 3, 1},
    {1, 3, 2},
    {1, 3, 3},
    {2, 1, 1},
    {2, 1, 2},
    {2, 1, 3},
    {2, 2, 1},
    {2, 2, 2},
    {2, 2, 3},
    {2, 3, 1},
    {2, 3, 2},
    {2, 3, 3},
    {3, 1, 1},
    {3, 1, 2},
    {3, 1, 3},
    {3, 2, 1},
    {3, 2, 2},
    {3, 3, 2},
    {3, 3, 1},
    {3, 3, 2},
    {3, 3, 3}
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hyperbolicConditionalMomentInversion::hyperbolicConditionalMomentInversion
(
    const dictionary& dict,
    const label nGeometricD
)
:
    nGeometricD_(nGeometricD),
    nDimensions_(3),
    nMoments_(nGeometricD_ == 2 ? 12 : 23),
    nNodes_(nGeometricD_ == 2 ? 9 : 27),
    support_("R"),
    moments_
    (
        nMoments_,
        nGeometricD_ == 2 ? twoDimMomentOrders_ : threeDimMomentOrders_
    ),
    abscissae_
    (
        nNodes_,
        nGeometricD_ == 2 ? twoDimNodeIndexes_ : threeDimNodeIndexes_
    ),
    weights_
    (
        nNodes_,
        nGeometricD_ == 2 ? twoDimNodeIndexes_ : threeDimNodeIndexes_
    ),
    univariateInverter_
    (
        new hyperbolicMomentInversion(dict.subDict("basicQuadrature"))
    ),
    etaMin_(dict.lookupOrDefault("etaMin", 1.0e-10)),
    qMax_(dict.lookupOrDefault("qMax", 30.0)),
    smallNegRealizability_
    (
        dict.lookupOrDefault
        (
            "smallNegRealizability",
            1.0e-6
        )
    ),
    varMin_(dict.lookupOrDefault("varMin", 1.0e-10))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hyperbolicConditionalMomentInversion
::~hyperbolicConditionalMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::hyperbolicConditionalMomentInversion::calcQ
(
    scalar q,
    scalar eta
)
{
    if (mag(q) > SMALL)
    {
        scalar slope = (eta - 3.0)/q;
        scalar sqrtDet = sqrt(8.0 + sqr(slope));

        if (q > 0.0)
        {
            return (slope + sqrtDet)/2.0;
        }
        else
        {
            return (slope - sqrtDet)/2.0;
        }
    }
    else
    {
        return 0.0;
    }
}

void Foam::hyperbolicConditionalMomentInversion::realizabilityUnivariateMoments
(
    scalar& c2,
    scalar& c3,
    scalar& c4
)
{
    if (c2 < 0.0)
    {
        c2 = 0.0;
        c3 = 0.0;
        c4 = 0.0;
    }

    if (c2*c4 < pow3(c2) + sqr(c3))
    {
        scalar q = c3/pow(c2, 3.0/2.0);
        scalar eta = c4/sqr(c2);
        q = calcQ(q, eta);
        eta = 1.0 + sqr(q);
        c3 = q*pow(c2, 3.0/2.0);
        c4 = eta*sqr(c2);
    }
}

void Foam::hyperbolicConditionalMomentInversion::invert1D
(
    const multivariateMomentSet& moments
)
{
    NotImplemented;
}

void Foam::hyperbolicConditionalMomentInversion::invert2D
(
    const multivariateMomentSet& moments,
    mappedList<scalar>& weights2D,
    mappedList<vector2D>& abscissae2D
)
{
    scalar m00 = moments(0, 0);

    if (m00 < SMALL)
    {
        forAll(weights2D, wi)
        {
            weights2D[wi] = m00;
            return;
        }
    };

    // Calculate normalized moments
    mappedScalarList scaledMoments(moments);

    forAll(scaledMoments, mi)
    {
        scaledMoments[mi] /= m00;
    }

    // Mean velocities and powers
    scalar meanU = scaledMoments(1, 0);
    scalar meanV = scaledMoments(0, 1);
    scalar sqrMeanU = sqr(meanU);
    scalar sqrMeanV = sqr(meanV);

    // Calculate central moments
    mappedScalarList centralMoments(scaledMoments);

    centralMoments(1, 0) = 0.0;
    centralMoments(0, 1) = 0.0;

    // NOTE: Sign changed due to -= operator
    centralMoments(2, 0) -= sqrMeanU;
    centralMoments(1, 1) -= meanU*meanV;
    centralMoments(0, 2) -= sqrMeanV;
    centralMoments(3, 0) -= 3.0*meanU*scaledMoments(2, 0) - 2.0*pow3(meanU);

    centralMoments(2, 1) -= meanV*scaledMoments(2, 0)
        + 2.0*meanU*scaledMoments(1, 1) - 2.0*sqrMeanU*meanV;

    centralMoments(1, 2) -= meanU*scaledMoments(0, 2)
        + 2.0*meanV*scaledMoments(1, 1) - 2.0*meanU*sqrMeanV;

    centralMoments(0, 3) -= 3.0*meanV*scaledMoments(0, 2) - 2.0*pow3(meanV);
    centralMoments(4, 0) -= 4.0*meanU*scaledMoments(3, 0)
        - 6.0*sqrMeanU*scaledMoments(2, 0) + 3.0*sqr(sqrMeanU);

    centralMoments(0, 4) -= 4.0*meanV*scaledMoments(0, 3)
        - 6.0*sqrMeanV*scaledMoments(0, 2) + 3.0*sqr(sqrMeanV);

    // One-dimensional inversion with realizability test
    univariateMomentSet mDir1
    (
        {
            1.0,
            0.0,
            centralMoments(2, 0),
            centralMoments(3, 0),
            centralMoments(4, 0)
        },
        "R"
    );

    // Find univariate quadrature in first direction
    univariateInverter_().invert(mDir1);

    // Store univariate quadrature in first direction
    scalarList wDir1(univariateInverter_().weights());
    scalarList absDir1(univariateInverter_().abscissae());

    // Default values of quadrature in direction 2
    forAll(weights2D, wi)
    {
        weights2D[wi] = 0.0;
        abscissae2D[wi] = vector2D::zero;
    }

    weights2D(1, 2) = 1.0;
    weights2D(3, 2) = 1.0;

    // Reconstruction settings
    bool linearVf = false;
    scalarList Vf(3, 0.0);

    if (centralMoments(2, 0) < varMin_)
    {
        // Resetting quadrature in the first direction
        wDir1 = 0.0;
        wDir1[1] = 1.0;

        univariateMomentSet mDir2
        (
            {
                1.0,
                0.0,
                centralMoments(0, 2),
                centralMoments(0, 3),
                centralMoments(0, 4)
            },
         "R"
        );

        //NOTE: Leave Vf elements null.

        univariateInverter_().invert(mDir1);

        for (label j = 1; j <= 3; j++)
        {
            weights2D(2, j) = univariateInverter_().weights()[j-1];
            abscissae2D(2, j).y() = univariateInverter_().abscissae()[j-1];
        }

    }
    else
    {
        scalar c20 = centralMoments(2, 0);
        scalar c02 = centralMoments(0, 2);
        scalar scaleFactor = sqrt(c20);
        scalarList scaledAbsDir1(absDir1);

        forAll(scaledAbsDir1, i)
        {
            scaledAbsDir1[i] /= scaleFactor;
        }

        scalar q = centralMoments(3, 0)/pow3(scaleFactor);
        scalar eta = centralMoments(4, 0)/sqr(c20);
        scalar scaledC11 = centralMoments(1, 1)/scaleFactor;
        scalar scaledC21 = centralMoments(2, 1)/c20;

        if (sqr(scaledC11) >= c02*(1.0 - varMin_))
        {
            scaledC11 = sign(scaledC11)*sqrt(c02);
            linearVf = true;
        }

        scalar realizability = eta - sqr(q) - 1.0;

        // Compute coefficients of function Vf
        scalarList VfCoeffs(2, 0.0);

        if (realizability > smallNegRealizability_ && !linearVf)
        {
            VfCoeffs[0] = (scaledC21 - q*scaledC11)/realizability;
            VfCoeffs[1] = ((eta - 1.0)*scaledC11 - q*scaledC21)/realizability;
        }
        else
        {
            VfCoeffs[1] = scaledC11;
        }

        forAll(Vf, vi)
        {
            Vf[vi] = VfCoeffs[1]*scaledAbsDir1[vi]
                + VfCoeffs[0]*(sqr(scaledAbsDir1[vi]) - 1.0);
        }

        // Find conditional variances
        scalarList bCoeffs(2, 0.0);

        forAll(Vf, vi)
        {
            bCoeffs[0] -= wDir1[vi]*sqr(Vf[vi]);
        }
        bCoeffs[0] += c02;

        if (bCoeffs[0] < 0.0)
        {
            if (!linearVf)
            {
                forAll(Vf, vi)
                {
                    Vf[vi] = scaledC11*scaledAbsDir1[vi];
                }

                bCoeffs[0] = 0.0;

                forAll(Vf, vi)
                {
                    bCoeffs[0] -= wDir1[vi]*sqr(Vf[vi]);
                }

                bCoeffs[0] += c02;
                bCoeffs[0] = max(0.0, bCoeffs[0]);

                linearVf = true;
            }
            else
            {
                bCoeffs[0] = 0.0;
            }
        }

        if (linearVf)
        {
            bCoeffs[1] = 0.0;
        }
        else
        {
            scalar scaledC12 = centralMoments(1, 2)/scaleFactor;

            forAll(Vf, vi)
            {
                bCoeffs[1] -= wDir1[vi]*scaledAbsDir1[vi]*sqr(Vf[vi]);
            }

            bCoeffs[1] += scaledC12;
        }

        scalarList mu2(3, 0.0);

        forAll(mu2, i)
        {
            mu2[i] = bCoeffs[0] + bCoeffs[1]*scaledAbsDir1[i];
        }

        //Find value and position of minimum value of mu
        label minMu2Index = findMin(mu2);
        scalar minMu2 = mu2[minMu2Index];

        if (minMu2 < 0.0)
        {
            bCoeffs[1] = -bCoeffs[0]/scaledAbsDir1[minMu2Index];

            forAll(mu2, i)
            {
                mu2[i] = bCoeffs[0] + bCoeffs[1]*scaledAbsDir1[i];
            }

            mu2[minMu2Index] = 0.0;
        }

        //Resetting q and eta
        q = 0.0;
        eta = 1.0;

        scalar sum1 = 0.0;
        forAll(wDir1, wi)
        {
            sum1 += wDir1[wi]*pow(mu2[wi], 3.0/2.0);
        }

        if (sum1 > SMALL)
        {
            scalar sum3 = 0.0;
            forAll(wDir1, wi)
            {
                sum3 -= wDir1[wi]*(pow3(Vf[wi]) + 3.0*Vf[wi]*mu2[wi]);
            }

            sum3 += centralMoments(0, 3);
            q = sum3/sum1;
        }

        scalar sum2 = 0.0;
        forAll(wDir1, wi)
        {
            sum2 += wDir1[wi]*sqr(mu2[wi]);
        }

        if (sum1 > SMALL)
        {
            scalar sum4a = 0.0;
            scalar sum4b = 0.0;

            forAll(wDir1, wi)
            {
                sum4a -= wDir1[wi]*pow4(Vf[wi])
                    + 6.0*wDir1[wi]*sqr(Vf[wi])*mu2[wi];

                sum4b -= 4.0*wDir1[wi]*Vf[wi]*pow(mu2[wi], 3.0/2.0);
            }

            sum4a += centralMoments(0, 4);
            eta = (sum4a + q*sum4b)/sum2;

            if (eta < sqr(q) + 1.0)
            {
                q = calcQ(q, eta);
                eta  = 1.0 + sqr(q);
            }
        }

        scalarList mu3(3, 0.0);
        scalarList mu4(3, 0.0);

        forAll(mu3, i)
        {
            mu3[i] = q*pow(mu2[i], 3.0/2.0);
            mu4[i] = eta*sqr(mu2[i]);
        }

        univariateMomentSet mMu1
        (
            {
                1.0,
                0.0,
                mu2[0],
                mu3[0],
                mu4[0]
            },
            "R"
        );

        univariateInverter_().invert(mMu1);

        for (label j = 1; j <= 3; j++)
        {
            weights2D(1, j) = univariateInverter_().weights()[j-1];
            abscissae2D(1, j).y() = univariateInverter_().abscissae()[j-1];
        }

        univariateMomentSet mMu2
        (
            {
                1.0,
                0.0,
                mu2[1],
                mu3[1],
                mu4[1]
            },
            "R"
        );

        univariateInverter_().invert(mMu2);

        for (label j = 1; j <= 3; j++)
        {
            weights2D(2, j) = univariateInverter_().weights()[j-1];
            abscissae2D(2, j).y() = univariateInverter_().abscissae()[j-1];
        }

        univariateMomentSet mMu3
        (
            {
                1.0,
                0.0,
                mu2[2],
                mu3[2],
                mu4[2]
            },
            "R"
        );

        univariateInverter_().invert(mMu3);

        for (label j = 1; j <= 3; j++)
        {
            weights2D(3, j) = univariateInverter_().weights()[j-1];
            abscissae2D(3, j).y() = univariateInverter_().abscissae()[j-1];
        }
    }

    // Compute multivariate quadrature
    for (label i = 1; i <= 3; i++)
    {
        for (label j = 1; j <= 3; j++)
        {
            weights2D(i, j) *= m00*wDir1[i];

            abscissae2D(i, j) =
                vector2D
                (
                    absDir1[i] + meanU,
                    Vf[i - 1]*abscissae2D(i, j).y() + meanV
                );
        }
    }
}

void Foam::hyperbolicConditionalMomentInversion::invert3D
(
    const multivariateMomentSet& moments
)
{
    scalar m000 = moments(0, 0, 0);

    if (m000 < SMALL)
    {
        forAll(weights_, wi)
        {
            weights_[wi] = m000;
            return;
        }
    };

    // Calculate normalized moments
    mappedScalarList scaledMoments(moments);

    forAll(scaledMoments, mi)
    {
        scaledMoments[mi] /= m000;
    }

    // Mean velocities and powers
    scalar meanU = scaledMoments(1, 0, 0);
    scalar meanV = scaledMoments(0, 1, 0);
    scalar meanW = scaledMoments(0, 0, 1);
    scalar sqrMeanU = sqr(meanU);
    scalar sqrMeanV = sqr(meanV);
    scalar sqrMeanW = sqr(meanW);

    // Calculate central moments
    mappedScalarList centralMoments(scaledMoments);

    centralMoments(1, 0, 0) = 0.0;
    centralMoments(0, 1, 0) = 0.0;
    centralMoments(0, 0, 1) = 0.0;

    // NOTE: Sign changed due to -= operator
    centralMoments(2, 0, 0) -= sqrMeanU;
    centralMoments(0, 2, 0) -= sqrMeanV;
    centralMoments(0, 0, 2) -= sqrMeanW;
    centralMoments(3, 0, 0) -= 3.0*meanU*scaledMoments(2, 0, 0)
        - 2.0*pow3(meanU);

    centralMoments(0, 3, 0) -= 3.0*meanV*scaledMoments(0, 2, 0)
        - 2.0*pow3(meanV);

    centralMoments(0, 0, 3) -= 3.0*meanW*scaledMoments(0, 0, 2)
        - 2.0*pow3(meanW);

    centralMoments(4, 0, 0) -= 4.0*meanU*scaledMoments(3, 0, 0)
        - 6.0*sqrMeanU*scaledMoments(2, 0, 0) + 3.0*sqr(meanU);

    centralMoments(0, 4, 0) -= 4.0*meanV*scaledMoments(0, 3, 0)
        - 6.0*sqrMeanV*scaledMoments(0, 2, 0) + 3.0*sqr(meanV);

    centralMoments(0, 0, 4) -= 4.0*meanW*scaledMoments(0, 0, 3)
        - 6.0*sqrMeanW*scaledMoments(0, 0, 3) + 3.0*sqr(meanW);

    if (m000 < SMALL)
    {
        // Use isotropic central moments
        centralMoments(1, 1, 0) = 0.0;
        centralMoments(1, 0, 1) = 0.0;
        centralMoments(0, 1, 1) = 0.0;
        centralMoments(2, 1, 0) = 0.0;
        centralMoments(2, 0, 1) = 0.0;
        centralMoments(1, 2, 0) = 0.0;
        centralMoments(1, 1, 1) = 0.0;
        centralMoments(1, 0, 2) = 0.0;
        centralMoments(0, 2, 1) = 0.0;
        centralMoments(0, 1, 2) = 0.0;
    }
    else
    {
        centralMoments(1, 1, 0) -= meanU*meanV;
        centralMoments(1, 0, 1) -= meanU*meanW;
        centralMoments(0, 1, 1) -= meanV*meanW;
        centralMoments(2, 1, 0) -= meanV*scaledMoments(2, 0, 0)
            + 2.0*meanU*scaledMoments(1, 1, 0) - 2.0*sqrMeanU*meanV;

        centralMoments(2, 0, 1) -= meanW*scaledMoments(2, 0, 0)
            + 2.0*meanU*scaledMoments(1, 0, 1) - 2.0*sqrMeanU*meanW;

        centralMoments(1, 2, 0) -= meanU*scaledMoments(1, 1, 0)
            - 2.0*meanU*sqrMeanV;

        centralMoments(1, 1, 1) -= meanV*scaledMoments(1, 0, 1)
            + meanW*scaledMoments(1, 1, 0) - 2.0*meanU*meanV*meanW;

        centralMoments(1, 0, 2) -= meanU*scaledMoments(0, 0, 2)
            + 2.0*meanW*scaledMoments(1, 0, 1) - 2.0*meanU*sqrMeanW;

        centralMoments(0, 2, 1) -= meanW*scaledMoments(0, 2, 0)
            + 2.0*meanV*scaledMoments(0, 1, 1) - 2.0*sqrMeanV*meanW;

        centralMoments(0, 1, 2) -= meanW*scaledMoments(0, 0, 2)
            + 2.0*meanW*scaledMoments(0, 1, 1) - 2.0*meanV*sqrMeanW;
    }

    // Check realizability of univariate moments
    realizabilityUnivariateMoments
    (
        centralMoments(2, 0, 0),
        centralMoments(3, 0, 0),
        centralMoments(4, 0, 0)
    );

    realizabilityUnivariateMoments
    (
        centralMoments(0, 2, 0),
        centralMoments(0, 3, 0),
        centralMoments(0, 4, 0)
    );

    realizabilityUnivariateMoments
    (
        centralMoments(0, 0, 2),
        centralMoments(0, 0, 3),
        centralMoments(0, 0, 4)
    );

    // One-dimensional inversion with realizability test
    univariateMomentSet mDir1
    (
        {
            1.0,
            0.0,
            centralMoments(2, 0, 0),
            centralMoments(3, 0, 0),
            centralMoments(4, 0, 0)
        },
        "R"
    );

    // Find univariate quadrature in first direction
    univariateInverter_().invert(mDir1);

    // Store univariate quadrature in first direction
    scalarList wDir1(univariateInverter_().weights());
    scalarList absDir1(univariateInverter_().abscissae());

    // Quadrature in second direction
    mappedList<scalar> w2DDir2(9, twoDimNodeIndexes_, 0.0);
    mappedList<vector2D> abs2Dir2(9, twoDimNodeIndexes_, vector2D::zero);

    // Default values
    w2DDir2(1, 2) = 1.0;
    w2DDir2(3, 2) = 1.0;

    // Default values of quadrature in direction 3
    forAll(weights_, wi)
    {
        weights_[wi] = 0.0;
        abscissae_[wi] = vector::zero;
    }

    weights_(1, 1, 2) = 1.0;
    weights_(1, 2, 2) = 1.0;
    weights_(1, 3, 2) = 1.0;
    weights_(2, 1, 2) = 1.0;
    weights_(2, 2, 2) = 1.0;
    weights_(2, 3, 2) = 1.0;
    weights_(3, 1, 2) = 1.0;
    weights_(3, 2, 2) = 1.0;
    weights_(3, 3, 2) = 1.0;
    weights_(3, 1, 2) = 1.0;

    scalarList VfCoeffs(3, 0.0);
    scalarSquareMatrix WfCoeffs(3, 0.0);

    if (centralMoments(2, 0, 0) < varMin_)
    {
        if (centralMoments(0, 2, 0) < varMin_)
        {
            univariateMomentSet mDir3
            (
                {
                    1.0,
                    0.0,
                    centralMoments(2, 0, 0),
                    centralMoments(3, 0, 0),
                    centralMoments(4, 0, 0)
                },
                "R"
            );

            // Find univariate quadrature in first direction
            univariateInverter_().invert(mDir3);

            // Store univariate quadrature in first direction
            scalarList wDir3(univariateInverter_().weights());
            scalarList absDir3(univariateInverter_().abscissae());

            wDir1[0] = 0.0;
            wDir1[1] = 1.0;
            wDir1[2] = 0.0;

            w2DDir2(2, 2) = 1.0;

            for (label wi = 1; wi <= 3; wi++)
            {
                absDir1[wi] = 0.0;
                weights_(2, 2, wi) = wDir3[wi - 1];
                abscissae_(2, 2, wi).z() = absDir3[wi - 1];
            }
        }
        else
        {
            multivariateMomentSet mDir23
            (
                {
                    1.0,
                    0.0,
                    0.0,
                    centralMoments(0, 2, 0),
                    centralMoments(0, 1, 1),
                    centralMoments(0, 0, 2),
                    centralMoments(0, 3, 0),
                    centralMoments(0, 2, 1),
                    centralMoments(0, 1, 2),
                    centralMoments(0, 0, 3),
                    centralMoments(0, 4, 0),
                    centralMoments(0, 0, 4)
                },
                twoDimMomentOrders_,
                "R"
            );

            invert2D(mDir23, w2DDir2, abs2Dir2);

            // TODO: Build quadrature here and return

        }
    }
}

void Foam::hyperbolicConditionalMomentInversion::invert
(
    const multivariateMomentSet& moments
)
{
    if (nGeometricD_ == 3)
    {
        //invert3D(moments);
    }
    else if (nGeometricD_ == 2)
    {
        //invert2D(moments);
    }
    else
    {
        invert1D(moments);
    }
}


// ************************************************************************* //
