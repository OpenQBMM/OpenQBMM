/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 Alberto Passalacqua
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

#include "CHyQMOMMomentInversion.H"
#include "mappedLists.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    defineTypeNameAndDebug(CHyQMOM, 0);
    addToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        CHyQMOM,
        dictionary
    );
}
}


const Foam::labelListList
Foam::multivariateMomentInversions::CHyQMOM::CHyQMOM::twoDimMomentOrders =
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
Foam::multivariateMomentInversions::CHyQMOM::CHyQMOM::threeDimMomentOrders =
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
Foam::multivariateMomentInversions::CHyQMOM::CHyQMOM::twoDimNodeIndexes =
{
    {0, 0},
    {0, 1},
    {0, 2},
    {1, 0},
    {1, 1},
    {1, 2},
    {2, 0},
    {2, 1},
    {2, 2}
};

const Foam::labelListList
Foam::multivariateMomentInversions::CHyQMOM::CHyQMOM::threeDimNodeIndexes =
{
    {0, 0, 0},
    {0, 0, 1},
    {0, 0, 2},
    {0, 1, 0},
    {0, 1, 1},
    {0, 1, 2},
    {0, 2, 0},
    {0, 2, 1},
    {0, 2, 2},
    {1, 0, 0},
    {1, 0, 1},
    {1, 0, 2},
    {1, 1, 0},
    {1, 1, 1},
    {1, 1, 2},
    {1, 2, 0},
    {1, 2, 1},
    {1, 2, 2},
    {2, 0, 0},
    {2, 0, 1},
    {2, 0, 2},
    {2, 1, 0},
    {2, 1, 1},
    {2, 1, 2},
    {2, 2, 0},
    {2, 2, 1},
    {2, 2, 2}
};


Foam::label Foam::multivariateMomentInversions::CHyQMOM::getNMoments
(
    const label nDims
)
{
    if (nDims == 1)
    {
        return 5;
    }
    else if (nDims == 2)
    {
        return 10;
    }
    else if (nDims == 3)
    {
        return 16;
    }
    return 0;
}


Foam::labelListList
Foam::multivariateMomentInversions::CHyQMOM::getMomentOrders
(
    const label nDims
)
{
    if (nDims == 1)
    {
        return {{0}, {1}, {2}, {3}, {4}};
    }
    else if (nDims == 2)
    {
        return twoDimMomentOrders;
    }
    else if (nDims == 3)
    {
        return threeDimMomentOrders;
    }

    return {{}};
}


Foam::label Foam::multivariateMomentInversions::CHyQMOM::getNNodes
(
    const label nDims
)
{
    if (nDims == 1)
    {
        return 3;
    }
    else if (nDims == 2)
    {
        return 9;
    }
    else if (nDims == 3)
    {
        return 27;
    }

    return 0;
}


Foam::labelListList
Foam::multivariateMomentInversions::CHyQMOM::getNodeIndexes
(
    const label nDims
)
{
    if (nDims == 1)
    {
        return {{0}, {1}, {2}};
    }
    else if (nDims == 2)
    {
        return twoDimNodeIndexes;
    }
    else if (nDims == 3)
    {
        return threeDimNodeIndexes;
    }

    return {{}};
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::CHyQMOM::CHyQMOM
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
:
    multivariateMomentInversion
    (
        dict,
        momentOrders,
        nodeIndexes,
        velocityIndexes
    ),
    univariateInverter_
    (
        new hyperbolicMomentInversion(dict)
    ),
    etaMin_(dict.lookupOrDefault<scalar>("etaMin", 1.0e-10)),
    qMax_(dict.lookupOrDefault<scalar>("qMax", 30.0)),
    smallNegRealizability_
    (
        dict.lookupOrDefault<scalar>
        (
            "smallNegRealizability",
            1.0e-6
        )
    ),
    varMin_(dict.lookupOrDefault<scalar>("varMin", 1.0e-10)),
    minCorrelation_(dict.lookupOrDefault<scalar>("minCorrelation", 1.0e-4))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::CHyQMOM::~CHyQMOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::multivariateMomentInversions::CHyQMOM::calcQ
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

void
Foam::multivariateMomentInversions::CHyQMOM::realizabilityUnivariateMoments
(
    scalar& c2,
    scalar& c3,
    scalar& c4
)
{
    if (c2 <= 0.0)
    {
        c2 = 0.0;
        c3 = 0.0;
        c4 = 0.0;

        return;
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

void Foam::multivariateMomentInversions::CHyQMOM::invert1D
(
    const multivariateMomentSet& moments,
    scalarList& weights1D,
    scalarList& abscissae1D
)
{
    scalar m0 = moments[0];
    label nWeights1D = weights1D.size();

    if (m0 < SMALL)
    {
        forAll(weights1D, wi)
        {
            weights1D[wi] = m0/scalar(nWeights1D);
            abscissae1D[wi] = 0.0;
        }

        return;
    };

    // Calculate normalized moments
    mappedScalarList scaledMoments(moments);

    forAll(scaledMoments, mi)
    {
        scaledMoments[mi] /= m0;
    }

    // Mean velocities and powers
    scalar meanU = scaledMoments(1);
    scalar sqrMeanU = sqr(meanU);

    // Calculate central moments
    mappedScalarList centralMoments(scaledMoments);

    centralMoments(1) = 0.0;

    // NOTE: Sign changed due to -= operator
    centralMoments(2) -= sqrMeanU;
    centralMoments(3) -= (3.0*meanU*scaledMoments(2) - 2.0*pow3(meanU));
    centralMoments(4) -= (4.0*meanU*scaledMoments(3)
        - 6.0*sqrMeanU*scaledMoments(2) + 3.0*sqr(sqrMeanU));

    // One-dimensional inversion with realizability test
    univariateMomentSet momentsToInvert
    (
        {
            scalar(1),
            scalar(0),
            centralMoments(2),
            centralMoments(3),
            centralMoments(4)
        },
        "R"
    );

    // Find univariate quadrature in first direction
    univariateInverter_().invert(momentsToInvert);

    // Store univariate quadrature in first direction
    forAll(weights1D, wi)
    {
        weights1D[wi] = m0*univariateInverter_().weights()[wi];
        abscissae1D[wi] = univariateInverter_().abscissae()[wi] + meanU;
    }
}

void Foam::multivariateMomentInversions::CHyQMOM::invert2D
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
    // Calculate normalized moments
    scalar s10 = moments(1, 0)/m00;
    scalar s01 = moments(0, 1)/m00;
    scalar s11 = moments(1, 1)/m00;
    scalar s20 = moments(2, 0)/m00;
    scalar s02 = moments(0, 2)/m00;
    scalar s30 = moments(3, 0)/m00;
    scalar s03 = moments(0, 3)/m00;
    scalar s40 = moments(4, 0)/m00;
    scalar s04 = moments(0, 4)/m00;

    // Mean velocities and powers
    scalar meanU = s10;
    scalar meanV = s01;
    scalar sqrMeanU = sqr(meanU);
    scalar sqrMeanV = sqr(meanV);

    // Calculate central moments
    scalar c20 = s20 - sqrMeanU;
    scalar c11 = s11 - meanU*meanV;
    scalar c02 = s02 - sqrMeanV;
    scalar c30 = s30;
    scalar c03 = s03;
    scalar c40 = s40;
    scalar c04 = s04;

    // NOTE: Sign changed due to -= operator
    c30 -= (3.0*meanU*s20 - 2.0*pow3(meanU));
    c03 -= (3.0*meanV*s02 - 2.0*pow3(meanV));
    c40 -= (4.0*meanU*s30 - 6.0*sqrMeanU*s20 + 3.0*sqr(sqrMeanU));
    c04 -= (4.0*meanV*s03 - 6.0*sqrMeanV*s02 + 3.0*sqr(sqrMeanV));

    realizabilityUnivariateMoments(c20, c30, c40);
    realizabilityUnivariateMoments(c02, c03, c04);
    // One-dimensional inversion with realizability test
    univariateMomentSet mDir1({scalar(1), scalar(0), c20, c30, c40}, "R");

    // Find univariate quadrature in first direction
    univariateInverter_->invert(mDir1);

    // Store univariate quadrature in first direction
    scalarList wDir1(univariateInverter_->weights());
    scalarList absDir1(univariateInverter_->abscissae());

    scalarList wDir2(3, Zero);
    scalarList absDir2(3, Zero);

    wDir2[1] = 1.0;

    scalarList Vf(3, Zero);

    if (c20 < varMin_)
    {
        wDir1 = 0.0;
        wDir1[1] = 1.0;

        univariateMomentSet mDir2({scalar(1), scalar(0), c02, c03, c04}, "R");

        //NOTE: Leave Vf elements null. AP
        univariateInverter_->invert(mDir2);

        forAll(wDir2, nodei)
        {
            wDir2[nodei] = univariateInverter_().weights()[nodei];
            absDir2[nodei] = univariateInverter_().abscissae()[nodei];
        }
    }
    else
    {
        scalar sqrtC20 = sqrt(c20);
        scalar c11s = c11/sqrtC20;

        if (sqr(c11s) > c02*(1.0 - 1e-10))
        {
            c11s = sign(c11s)*sqrt(c02);
        }

        forAll(Vf, vi)
        {
            Vf[vi] = c11s*absDir1[vi]/sqrtC20;
        }

        // Compute conditional variance
        scalar sumVars = 0.0;

        forAll(Vf, vi)
        {
            sumVars += wDir1[vi]*sqr(Vf[vi]);
        }

        scalar mu2Avg = max(c02 - sumVars, scalar(0));

        scalarList mu(5, Zero);
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

            scalar q = (c03 - sumWVf3)/sqrt(pow3(mu[2]));
            scalar eta = (c04 - sumWVf4 - 6.0*sumVars*mu[2])/sqr(mu[2]);

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
                scalar(1),
                scalar(0),
                mu[2],
                mu[3],
                mu[4]
            },
            "R"
        );

        univariateInverter_->invert(mMu);

        forAll(wDir2, nodei)
        {
            wDir2[nodei] = univariateInverter_().weights()[nodei];
            absDir2[nodei] = univariateInverter_().abscissae()[nodei];
        }

    }

    // Compute multivariate quadrature
    for (label i = 0; i < 3; i++)
    {
        for (label j = 0; j < 3; j++)
        {
            weights2D(i, j) = m00*wDir1[i]*wDir2[j];

            abscissae2D(i, j) =
                vector2D
                (
                    absDir1[i] + meanU,
                    Vf[i] + absDir2[j] + meanV
                );
        }
    }
}

void Foam::multivariateMomentInversions::CHyQMOM::invert3D
(
    const multivariateMomentSet& moments
)
{
    scalar m000 = moments(0, 0, 0);

    if (m000 < SMALL)
    {
        weights_(2,2,2) = m000;
        return;
    };

    // Calculate normalized moments
    scalar s100 = moments(1)/m000;
    scalar s010 = moments(0, 1)/m000;
    scalar s001 = moments(0, 0, 1)/m000;
    scalar s200 = moments(2)/m000;
    scalar s110 = moments(1, 1, 0)/m000;
    scalar s101 = moments(1, 0, 1)/m000;
    scalar s020 = moments(0, 2)/m000;
    scalar s011 = moments(0, 1, 1)/m000;
    scalar s002 = moments(0, 0, 2)/m000;
    scalar s300 = moments(3)/m000;
    scalar s030 = moments(0, 3)/m000;
    scalar s003 = moments(0, 0, 3)/m000;
    scalar s400 = moments(4)/m000;
    scalar s040 = moments(0, 4)/m000;
    scalar s004 = moments(0, 0, 4)/m000;

    // Calculate central moments
    scalar meanU = s100;
    scalar meanV = s010;
    scalar meanW = s001;
    scalar sqrMeanU = sqr(meanU);
    scalar sqrMeanV = sqr(meanV);
    scalar sqrMeanW = sqr(meanW);

    scalar c200 = s200;
    scalar c110 = s110;
    scalar c101 = s101;
    scalar c011 = s011;
    scalar c020 = s020;
    scalar c002 = s002;
    scalar c300 = s300;
    scalar c030 = s030;
    scalar c003 = s003;
    scalar c400 = s400;
    scalar c040 = s040;
    scalar c004 = s004;

    // NOTE: Sign changed due to -= operator
    c200 -= sqrMeanU;
    c110 -= meanU*meanV;
    c101 -= meanU*meanW;
    c020 -= sqrMeanV;
    c011 -= meanV*meanW;
    c002 -= sqrMeanW;

    c300 -= (3.0*meanU*s200 - 2.0*pow3(meanU));
    c030 -= (3.0*meanV*s020 - 2.0*pow3(meanV));
    c003 -= (3.0*meanW*s002 - 2.0*pow3(meanW));

    c400 -= (4.0*meanU*s300 - 6.0*sqrMeanU*s200 + 3.0*sqr(sqrMeanU));
    c040 -= (4.0*meanV*s030 - 6.0*sqrMeanV*s020 + 3.0*sqr(sqrMeanV));
    c004 -= (4.0*meanW*s003 - 6.0*sqrMeanW*s002 + 3.0*sqr(sqrMeanW));

    if (c200 <= 0.0)
    {
        c200 = 0.0;
        c300 = 0.0;
        c400 = 0.0;
    }

    if (c200*c400 < pow3(c200) + sqr(c300))
    {
        scalar q = c300/sqrt(pow3(c200));
        scalar eta = c400/sqr(c200);

        if (mag(q) > SMALL)
        {
            q = calcQ(q, eta);
            eta  = 1.0 + sqr(q);
        }

        c300 = q*sqrt(pow3(c200));
        c400 = eta*sqr(c200);
    }

    if (c020 <= 0.0)
    {
        c020 = 0.0;
        c030 = 0.0;
        c040 = 0.0;
    }

    if (c020*c040 < pow3(c020) + sqr(c030))
    {
        scalar q = c030/sqrt(pow3(c020));
        scalar eta = c040/sqr(c020);

        if (mag(q) > SMALL)
        {
            q = calcQ(q, eta);
            eta  = 1.0 + sqr(q);
        }

        c030 = q*sqrt(pow3(c020));
        c040 = eta*sqr(c020);
    }

    if (c002 <= 0.0)
    {
        c002 = 0.0;
        c003 = 0.0;
        c004 = 0.0;
    }

    if (c002*c004 < pow3(c002) + sqr(c003))
    {
        scalar q = c003/sqrt(pow3(c002));
        scalar eta = c004/sqr(c002);

        if (mag(q) > SMALL)
        {
            q = calcQ(q, eta);
            eta  = 1.0 + sqr(q);
        }

        c003 = q*sqrt(pow3(c002));
        c004 = eta*sqr(c002);
    }

    // Invert first direction
    univariateMomentSet mDir1({scalar(1), scalar(0), c200, c300, c400}, "R");

    // Find univariate quadrature in first direction
    univariateInverter_().invert(mDir1);

    // Store univariate quadrature in first direction
    scalarList wDir1(univariateInverter_().weights());
    scalarList absDir1(univariateInverter_().abscissae());

    // X direction is degenerate
    if (c200 < varMin_)
    {
        // X and y directions are degenerate
        if (c020 < varMin_)
        {
            univariateMomentSet mDir3({scalar(1), scalar(0), c002, c003, c004}, "R");

            // Find univariate quadrature in first direction
            univariateInverter_().invert(mDir3);

            // Store univariate quadrature in first direction
            scalarList wDir3(univariateInverter_().weights());
            scalarList absDir3(univariateInverter_().abscissae());

            for(label i = 0; i < 3; i++)
            {
                weights_(1, 1, i) = m000*wDir3[i];
                velocityAbscissae_(1, 1, i) =
                    vector
                    (
                        meanU,
                        meanV,
                        absDir3[i] + meanW
                    );
            }

            return;
        }
        // Only x direction is degenerate
        else
        {
            multivariateMomentSet mDir23
            (
                {
                    scalar(1),
                    scalar(0),
                    scalar(0),
                    c020,
                    c011,
                    c002,
                    c030,
                    c003,
                    c040,
                    c004
                },
                twoDimMomentOrders,
                "R"
            );

            mappedList<scalar> wDir23(9, twoDimNodeIndexes, Zero);
            mappedList<vector2D> absDir23
            (
                9,
                twoDimNodeIndexes,
                vector2D::zero
            );

            invert2D(mDir23, wDir23, absDir23);

            for (label j = 0; j < 3; j++)
            {
                for (label k = 0; k < 3; k++)
                {
                    weights_(1, j, k) = m000*wDir23(j, k);
                    velocityAbscissae_(1, j, k) =
                        vector
                        (
                            meanU,
                            absDir23(j, k).x() + meanV,
                            absDir23(j, k).y() + meanW
                        );
                }
            }

            return;
        }
    }
    // Y direction is degenerate
    else if (c020 < varMin_)
    {
        multivariateMomentSet mDir13
        (
            {
                scalar(1),
                scalar(0),
                scalar(0),
                c200,
                c101,
                c002,
                c300,
                c003,
                c400,
                c004
            },
            twoDimMomentOrders,
           "R"
        );

        mappedList<scalar> wDir13(9, twoDimNodeIndexes, Zero);
        mappedList<vector2D> absDir13
        (
            9,
            twoDimNodeIndexes,
            Zero
        );

        invert2D(mDir13, wDir13, absDir13);

        for (label i = 0; i < 3; i++)
        {
            for (label k = 0; k < 3; k++)
            {
                weights_(i, 1, k) = m000*wDir13(i, k);
                velocityAbscissae_(i, 1, k) =
                    vector
                    (
                        absDir13(i, k).x() + meanU,
                        meanV,
                        absDir13(i, k).y() + meanW
                    );
            }
        }

        return;
    }
    // X and y directions are non-degenerate
    else
    {
        multivariateMomentSet mDir12
        (
            {
                scalar(1),
                scalar(0),
                scalar(0),
                c200,
                c110,
                c020,
                c300,
                c030,
                c400,
                c040
            },
            twoDimMomentOrders,
           "R"
        );

        mappedList<scalar> wDir12(9, twoDimNodeIndexes, Zero);
        mappedList<vector2D> abscissaeDir12
        (
            9,
            twoDimNodeIndexes,
            Zero
        );

        invert2D(mDir12, wDir12, abscissaeDir12);

        // Z direction is degenerate
        if (c002 < varMin_)
        {
            for (label i = 0; i < 3; i++)
            {
                for (label j = 0; j < 3; j++)
                {
                    weights_(i, j, 1) = m000*wDir12(i, j);
                    velocityAbscissae_(i, j, 1).x() =
                        abscissaeDir12(i, j).x() + meanU;
                    velocityAbscissae_(i, j, 1).y() =
                        abscissaeDir12(i, j).y() + meanV;
                }
            }

            return;
        }
        // All directions are non-degenerate
        else
        {
            // Scale weights in directions 12
            scalar sumWeights1 = 0.0;
            scalar sumWeights2 = 0.0;
            scalar sumWeights3 = 0.0;

            for (label i = 0; i < 3; i++)
            {
                sumWeights1 += wDir12(0, i);
                sumWeights2 += wDir12(1, i);
                sumWeights3 += wDir12(2, i);
            }

            for (label i = 0; i < 3; i++)
            {
                wDir12(0, i) /= sumWeights1;
                wDir12(1, i) /= sumWeights2;
                wDir12(2, i) /= sumWeights3;
            }

            // Compute Vf reconstruction
            scalarList Vf(3, Zero);
            for (label i = 0; i < 3; i++)
            {
                Vf[0] += wDir12(0, i)*abscissaeDir12(0, i).y();
                Vf[1] += wDir12(1, i)*abscissaeDir12(1, i).y();
                Vf[2] += wDir12(2, i)*abscissaeDir12(2, i).y();
            }

            mappedList<scalar> absDir12(9, twoDimNodeIndexes, Zero);
            for (label i = 0; i < 3; i++)
            {
                for (label j = 0; j < 3; j++)
                {
                    absDir12(i, j) = abscissaeDir12(i, j).y() - Vf[i];
                }
            }

            scalar sqrtC200 = sqrt(c200);
            scalar sqrtC020 = sqrt(c020);
            scalar sqrtC002 = sqrt(c002);

            scalarSquareMatrix W1(3, 0.0);
            scalarSquareMatrix W2(3, 0.0);
            scalarSquareMatrix Vp(3, 0.0);

            for (label i = 0; i < 3; i++)
            {
                W1(i, i) = wDir1[i];
                for (label j = 0; j < 3; j++)
                {
                    W2(i, j) = wDir12(i, j);
                    Vp(i, j) = absDir12(i, j);
                }
            }

            scalarSquareMatrix RAB = W1*W2;

            //scalarSquareMatrix Vc0(3, 1.0);
            scalarSquareMatrix Vc1(3, 0.0);
            scalarSquareMatrix Vc2(3, 0.0);
            scalarSquareMatrix C01(3, 0.0);

            scalar A1 = 0.0;
            scalar A2 = 0.0;

            for (label i = 0; i < 3; i++)
            {
                for (label j = 0; j < 3; j++)
                {
                    Vc1(i, j) = absDir1[i]/sqrtC200;
                    Vc2(i, j) = Vp(i, j)/sqrtC020;

                    scalar VAB = (Vp(i,j) + Vf[i])/sqrtC020;
                    C01(i, j) += RAB(i,j)*VAB;

                    A1 += C01(i, j)*Vc1(i, j);
                    A2 += C01(i, j)*Vc2(i, j);
                }
            }

            scalar c101s = c101/sqrtC200;
            scalar c011s = c011/sqrtC020;

            if (sqr(c101s) >= c002*(1.0 - SMALL))
            {
                c101s = sign(c101s)*sqrtC002;
            }
            else if (sqr(c011s) >= c002*(1.0 - SMALL))
            {
                scalar c110s = c110/(sqrtC200*sqrtC020);
                c011s = sign(c011s)*sqrtC002;
                c101s = c110s*c011s;
            }

            //scalar b0 = 0;
            scalar b1 = c101s;
            scalar b2 = 0;

            if (A2 > minCorrelation_)
            {
                b2 = (c011s - A1*b1)/A2;
            }

            scalar sumVarDir3 = 0.0;
            scalarSquareMatrix Wf(3, 0.0);
            for (label i = 0; i < 3; i++)
            {
                for (label j = 0; j < 3; j++)
                {
                    Wf(i, j) = b1*Vc1(i, j) + b2*Vc2(i, j);// + b0*Vc0(i, j)
                    sumVarDir3 += RAB(i, j)*sqr(Wf(i, j));
                }
            }

            scalarList mu(5, Zero);
            mu[2] = max(c002 - sumVarDir3, scalar(0));

            scalar q = 0;
            scalar eta = 1;
            if (mu[2] > etaMin_)
            {
                scalar sum1 = mu[2]*sqrt(mu[2]);
                scalar sum2 = sqr(mu[2]);
                scalar sum3 = 0.0;
                scalar sum4 = 0.0;
                for (label i = 0; i < 3; i++)
                {
                    for (label j = 0; j < 3; j++)
                    {
                        sum3 += RAB(i, j)*pow3(Wf(i, j));
                        sum4 += RAB(i, j)*pow4(Wf(i, j));
                    }
                }
                sum4 += + 6.0*sumVarDir3*mu[2];

                q = (c003 - sum3)/sum1;
                eta = (c004 - sum4)/sum2;

                if (eta < sqr(q) + 1)
                {
                    q = calcQ(q, eta);
                    eta = sqr(q) + 1.0;
                }
            }
            mu[3] = q*mu[2]*sqrt(mu[2]);
            mu[4] = eta*sqr(mu[2]);

            // Invert final direction
            univariateMomentSet mDir3({scalar(1), scalar(0), mu[2], mu[3], mu[4]}, "R");

            // Find univariate quadrature in final direction
            univariateInverter_().invert(mDir3);

            scalarList wDir3(univariateInverter_().weights());
            scalarList absDir3(univariateInverter_().abscissae());

            for (label i = 0; i < 3; i++)
            {
                for (label j = 0; j < 3; j++)
                {
                    for (label k = 0; k < 3; k++)
                    {
                        weights_(i, j, k) =
                            m000*wDir1[i]*wDir12(i, j)*wDir3[k];

                        velocityAbscissae_(i, j, k) =
                            vector
                            (
                                absDir1[i] + meanU,
                                Vf[i] + absDir12(i, j) + meanV,
                                Wf(i, j) + absDir3[k] + meanW
                            );
                    }
                }
            }

        }
    }
}

bool Foam::multivariateMomentInversions::CHyQMOM::invert
(
    const multivariateMomentSet& moments
)
{
    reset();

    if (nvelocityDimensions_ == 3)
    {
        invert3D(moments);
    }
    else if (nvelocityDimensions_ == 2)
    {
        mappedScalarList w
        (
            getNNodes(2),
            twoDimNodeIndexes
        );
        mappedList<vector2D> u
        (
            getNNodes(2),
            twoDimNodeIndexes
        );

        invert2D(moments, w, u);

        forAll(u, nodei)
        {
            weights_[nodei] = w[nodei];
            velocityAbscissae_[nodei] =
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
        scalarList w(getNNodes(1), Zero);
        scalarList u(getNNodes(1), Zero);

        invert1D(moments, w, u);

        forAll(w, nodei)
        {
            weights_[nodei] = w[nodei];
            velocityAbscissae_[nodei] = vector(u[nodei], 0.0, 0.0);
        }
    }

    return true;
}


// ************************************************************************* //
