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
    {0, 3, 0},
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
    {3, 2, 3},
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
    nMoments_(nGeometricD_ == 2 ? 10 : 16),
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
    varMin_(dict.lookupOrDefault("varMin", 1.0e-10)),
    minCorrelation_(dict.lookupOrDefault("minCorrelation", 1.0e-4))
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
    label nWeights2D = weights2D.size();

    if (m00 < SMALL)
    {
        forAll(weights2D, wi)
        {
            weights2D[wi] = m00/scalar(nWeights2D);
            abscissae2D[wi] = vector2D::zero;
        }
        return;
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
    centralMoments(3, 0) -= (3.0*meanU*scaledMoments(2, 0) - 2.0*pow3(meanU));
    centralMoments(0, 3) -= (3.0*meanV*scaledMoments(0, 2) - 2.0*pow3(meanV));
    centralMoments(4, 0) -= (4.0*meanU*scaledMoments(3, 0)
        - 6.0*sqrMeanU*scaledMoments(2, 0) + 3.0*sqr(sqrMeanU));

    centralMoments(0, 4) -= (4.0*meanV*scaledMoments(0, 3)
        - 6.0*sqrMeanV*scaledMoments(0, 2) + 3.0*sqr(sqrMeanV));

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

    scalarList wDir2(3, 0.0);
    scalarList absDir2(3, 0.0);

    // Default values of quadrature in direction 2
    forAll(wDir2, nodei)
    {
        wDir2[nodei] = 0.0;
        absDir2[nodei] = 0.0;
    }

    wDir2[1] = 1.0;

    // Reconstruction settings
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

        //NOTE: Leave Vf elements null. AP

        univariateInverter_().invert(mDir2);

        forAll(wDir2, nodei)
        {
            wDir2[nodei] = univariateInverter_().weights()[nodei];
            absDir2[nodei] = univariateInverter_().abscissae()[nodei];
        }
    }
    else
    {
        forAll(Vf, vi)
        {
            Vf[vi] = centralMoments(1, 1)*absDir1[vi]/centralMoments(2, 0);
        }

        // Compute conditional variance
        scalar sumVars = 0.0;

        forAll(Vf, vi)
        {
            sumVars += wDir1[vi]*sqr(Vf[vi]);
        }

        scalar mu2Avg = max(centralMoments(0, 2) - sumVars, 0.0);

        scalarList mu(5, 0.0);
        mu[0] = 1.0;
        mu[2] = mu2Avg;
        mu[4] = sqr(mu2Avg);

        if (mu[2] > varMin_)
        {
            scalar sumWVf3 = 0.0;
            scalar sumWVf4 = 0.0;

            forAll(Vf, vi)
            {
                sumWVf3 += wDir1[vi]*pow3(Vf[vi]);
                sumWVf4 += wDir1[vi]*pow4(Vf[vi]);
            }

            scalar q = (centralMoments(0, 3) - sumWVf3)/sqrt(pow3(mu[2]));

            scalar eta =
                (centralMoments(0, 4) - sumWVf4 - 6.0*sumVars*mu[2])/sqr(mu[2]);

            if (eta < sqr(q) + 1.0)
            {
                q = calcQ(q, eta);
                eta  = 1.0 + sqr(q);
            }

            mu[3] = q*sqrt(pow3(mu[2]));
            mu[4] = eta*sqr(mu[2]);
        }

        univariateMomentSet mMu
        (
            {
                1.0,
                0.0,
                mu[2],
                mu[3],
                mu[4]
            },
            "R"
        );

        univariateInverter_().invert(mMu);

        forAll(wDir2, nodei)
        {
            wDir2[nodei] = univariateInverter_().weights()[nodei];
            absDir2[nodei] = univariateInverter_().abscissae()[nodei];
        }

    }

    // Compute multivariate quadrature
    for (label i = 1; i <= 3; i++)
    {
        for (label j = 1; j <= 3; j++)
        {
            weights2D(i, j) = m00*wDir1[i - 1]*wDir2[j - 1];

            abscissae2D(i, j) =
                vector2D
                (
                    absDir1[i - 1] + meanU,
                    Vf[i - 1] + absDir2[j-1] + meanV
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
    }
    else
    {
        centralMoments(1, 1, 0) -= meanU*meanV;
        centralMoments(1, 0, 1) -= meanU*meanW;
        centralMoments(0, 1, 1) -= meanV*meanW;
    }

    if (centralMoments(2, 0, 0) <= 0.0)
    {
        centralMoments(2, 0, 0) = 0.0;
        centralMoments(3, 0, 0) = 0.0;
        centralMoments(4, 0, 0) = 0.0;
    }

    scalar c200 = centralMoments(2, 0, 0);

    if
    (
        c200*centralMoments(4, 0, 0)
      < pow3(c200) + sqr(centralMoments(3, 0, 0))
    )
    {
        scalar q = centralMoments(3, 0, 0)/sqrt(pow3(c200));
        scalar eta = centralMoments(4, 0, 0)/sqr(c200);

        if (mag(q) > SMALL)
        {
            q = calcQ(q, eta);
            eta  = 1.0 + sqr(q);
        }

        centralMoments(3, 0, 0) = q*sqrt(pow3(c200));
        centralMoments(4, 0, 0) = eta*sqr(c200);
    }

    if (centralMoments(0, 2, 0) <= 0.0)
    {
        centralMoments(0, 2, 0) = 0.0;
        centralMoments(0, 3, 0) = 0.0;
        centralMoments(0, 4, 0) = 0.0;
    }

    scalar c020 = centralMoments(0, 2, 0);

    if
    (
        c020*centralMoments(0, 4, 0)
      < pow3(c020) + sqr(centralMoments(0, 3, 0))
    )
    {
        scalar q = centralMoments(0, 3, 0)/sqrt(pow3(c020));
        scalar eta = centralMoments(0, 4, 0)/sqr(c020);

        if (mag(q) > SMALL)
        {
            q = calcQ(q, eta);
            eta  = 1.0 + sqr(q);
        }

        centralMoments(0, 3, 0) = q*sqrt(pow3(c020));
        centralMoments(0, 4, 0) = eta*sqr(c020);
    }

    if (centralMoments(0, 2, 0) <= 0.0)
    {
        centralMoments(0, 2, 0) = 0.0;
        centralMoments(0, 3, 0) = 0.0;
        centralMoments(0, 4, 0) = 0.0;
    }

    scalar c002 = centralMoments(0, 0, 2);

    if
    (
        c002*centralMoments(0, 0, 4)
      < pow3(c002) + sqr(centralMoments(0, 0, 3))
    )
    {
        scalar q = centralMoments(0, 0, 3)/sqrt(pow3(c002));
        scalar eta = centralMoments(0, 0, 4)/sqr(c002);

        if (mag(q) > SMALL)
        {
            q = calcQ(q, eta);
            eta  = 1.0 + sqr(q);
        }

        centralMoments(0, 0, 3) = q*sqrt(pow3(c002));
        centralMoments(0, 0, 4) = eta*sqr(c002);
    }

    scalarList Vf(3, 0.0);
    scalarSquareMatrix Wf(3, 0.0);

    if (centralMoments(2, 0, 0) < varMin_)
    {
        if (centralMoments(0, 2, 0) < varMin_)
        {
            univariateMomentSet mDir3
            (
                {
                    1.0,
                    0.0,
                    centralMoments(0, 0, 2),
                    centralMoments(0, 0, 3),
                    centralMoments(0, 0, 4)
                },
                "R"
            );

            // Find univariate quadrature in first direction
            univariateInverter_().invert(mDir3);

            // Store univariate quadrature in first direction
            scalarList weightsDir3(univariateInverter_().weights());
            scalarList abscissaeDir3(univariateInverter_().abscissae());

            for(label i = 1; i <= 3; i++)
            {
                weights_(2, 2, i) = m000*weightsDir3[i - 1];
                abscissae_(2, 2, i).z() = abscissaeDir3[i - 1] + meanW;
            }

            return;
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
                    centralMoments(0, 0, 3),
                    centralMoments(0, 4, 0),
                    centralMoments(0, 0, 4)
                },
                twoDimMomentOrders_,
                "R"
            );

            mappedList<scalar> weightsDir23(9, twoDimNodeIndexes_, 0.0);
            mappedList<vector2D> abscissaeDir23
            (
                9,
                twoDimNodeIndexes_,
                vector2D::zero
            );

            invert2D(mDir23, weightsDir23, abscissaeDir23);

            for (label j = 1; j <= 3; j++)
            {
                for (label k = 1; k <= 3; k++)
                {
                    weights_(2, j, k) = m000*weightsDir23(j, k);
                    abscissae_(2, j, k).y() =
                        abscissaeDir23(j, k).x() + meanV;
                    abscissae_(2, j, k).z() =
                        abscissaeDir23(j, k).y() + meanW;
                }
            }

            return;
        }
    }
    else if (centralMoments(0, 2, 0) < varMin_)
    {
        multivariateMomentSet mDir13
        (
            {
                1.0,
                0.0,
                0.0,
                centralMoments(2, 0, 0),
                centralMoments(1, 0, 1),
                centralMoments(0, 0, 2),
                centralMoments(3, 0, 0),
                centralMoments(0, 0, 3),
                centralMoments(4, 0, 0),
                centralMoments(0, 0, 4)
            },
            twoDimMomentOrders_,
           "R"
        );

        mappedList<scalar> weightsDir13(9, twoDimNodeIndexes_, 0.0);
        mappedList<vector2D> abscissaeDir13
        (
            9,
            twoDimNodeIndexes_,
            vector2D::zero
        );

        invert2D(mDir13, weightsDir13, abscissaeDir13);

        for (label i = 1; i <= 3; i++)
        {
            for (label k = 1; k <= 3; k++)
            {
                weights_(i, 2, k) = m000*weightsDir13(i, k);
                abscissae_(i, 2, k).x() = abscissaeDir13(i, k).x() + meanU;
                abscissae_(i, 2, k).z() = abscissaeDir13(i, k).y() + meanW;
            }
        }

        return;
    }
    else
    {
        multivariateMomentSet mDir12
        (
            {
                1.0,
                0.0,
                0.0,
                centralMoments(2, 0, 0),
                centralMoments(1, 1, 0),
                centralMoments(0, 2, 0),
                centralMoments(3, 0, 0),
                centralMoments(0, 3, 0),
                centralMoments(4, 0, 0),
                centralMoments(0, 4, 0)
            },
            twoDimMomentOrders_,
           "R"
        );

        mappedList<scalar> weightsDir12(9, twoDimNodeIndexes_, 0.0);
        mappedList<vector2D> abscissaeDir12
        (
            9,
            twoDimNodeIndexes_,
            vector2D::zero
        );

        invert2D(mDir12, weightsDir12, abscissaeDir12);

        if (centralMoments(0, 0, 2) < varMin_)  // Degenerate third direction
        {
            for (label i = 1; i <= 3; i++)
            {
                for (label j = 1; j <= 3; j++)
                {
                    weights_(i, j, 2) = m000*weightsDir12(i, j);
                    abscissae_(i, j, 2).x() = abscissaeDir12(i, j).x() + meanU;
                    abscissae_(i, j, 2).y() = abscissaeDir12(i, j).y() + meanV;
                }
            }

            return;
        }
        else  // All directions are non-degenerate
        {
            Info<<"here"<<endl;
            // One-dimensional inversion
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
            scalarList weightsDir1(univariateInverter_().weights());
            scalarList abscissaeDir1(univariateInverter_().abscissae());

            // Scale weights in directions 12
            scalar sumWeights1 = 0.0;
            scalar sumWeights2 = 0.0;
            scalar sumWeights3 = 0.0;

            for (label i = 1; i <= 3; i++)
            {
                sumWeights1 += weightsDir12(1, i);
                sumWeights2 += weightsDir12(2, i);
                sumWeights3 += weightsDir12(3, i);
            }

            for (label i = 1; i <= 3; i++)
            {
                weightsDir12(1, i) /= sumWeights1;
                weightsDir12(2, i) /= sumWeights2;
                weightsDir12(3, i) /= sumWeights3;
            }

            // Compute Vf reconstruction
            forAll(Vf, vfi)
            {
                Vf[0] +=
                    weightsDir12(1, vfi + 1)*abscissaeDir12(1, vfi + 1).y();

                Vf[1] +=
                    weightsDir12(2, vfi + 1)*abscissaeDir12(2, vfi + 1).y();

                Vf[2] +=
                    weightsDir12(3, vfi + 1)*abscissaeDir12(3, vfi + 1).y();
            }
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
        invert3D(moments);
    }
    else if (nGeometricD_ == 2)
    {
        mappedList<scalar> w
        (
            nNodes_,
            twoDimNodeIndexes_
        );
        mappedList<vector2D> u
        (
            nNodes_,
            twoDimNodeIndexes_
        );

        invert2D(moments, w, u);

        forAll(u, nodei)
        {
            weights_[nodei] = w[nodei];
            abscissae_[nodei] =
                vector
                (
                    u[nodei].x(),
                    u[nodei].y(),
                    0.0
                );
        }
    }
    else
    {
        invert1D(moments);
    }
}


// ************************************************************************* //
