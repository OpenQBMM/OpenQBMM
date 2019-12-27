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

#include "CHyQMOMPlusMomentInversion.H"
#include "QRMatrix.H"
#include "SVD.H"
#include "mappedLists.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    defineTypeNameAndDebug(CHyQMOMPlus, 0);
    addToRunTimeSelectionTable
    (
        multivariateMomentInversion,
        CHyQMOMPlus,
        dictionary
    );
}
}


const Foam::labelListList
Foam::multivariateMomentInversions::CHyQMOMPlus::CHyQMOMPlus::twoDimMomentOrders =
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
Foam::multivariateMomentInversions::CHyQMOMPlus::CHyQMOMPlus::threeDimMomentOrders =
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
Foam::multivariateMomentInversions::CHyQMOMPlus::CHyQMOMPlus::twoDimNodeIndexes =
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
Foam::multivariateMomentInversions::CHyQMOMPlus::CHyQMOMPlus::threeDimNodeIndexes =
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


Foam::label Foam::multivariateMomentInversions::CHyQMOMPlus::getNMoments
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
        return 12;
    }
    else if (nDims == 3)
    {
        return 23;
    }
    return 0;
}


Foam::labelListList
Foam::multivariateMomentInversions::CHyQMOMPlus::getMomentOrders
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


Foam::label Foam::multivariateMomentInversions::CHyQMOMPlus::getNNodes
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
Foam::multivariateMomentInversions::CHyQMOMPlus::getNodeIndexes
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

Foam::multivariateMomentInversions::CHyQMOMPlus::CHyQMOMPlus
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
    qMax_(dict.lookupOrDefault<scalar>("qMax", 100.0)),
    smallNegRealizability_
    (
        dict.lookupOrDefault<scalar>
        (
            "smallNegRealizability",
            1.0e-6
        )
    ),
    varMin_(dict.lookupOrDefault<scalar>("varMin", 1.0e-10)),
    minCorrelation_(dict.lookupOrDefault<scalar>("minCorrelation", 1.0e-6))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversions::CHyQMOMPlus::~CHyQMOMPlus()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::multivariateMomentInversions::CHyQMOMPlus::calcQ
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
Foam::multivariateMomentInversions::CHyQMOMPlus::realizabilityUnivariateMoments
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

void Foam::multivariateMomentInversions::CHyQMOMPlus::invert1D
(
    const multivariateMomentSet& moments,
    scalarList& weights1D,
    scalarList& abscissae1D
)
{
    // One-dimensional inversion with realizability test
    univariateMomentSet momentsToInvert
    (
        {
            moments(0),
            moments(1),
            moments(2),
            moments(3),
            moments(4)
        },
        "R"
    );

    // Find univariate quadrature in first direction
    univariateInverter_().invert(momentsToInvert);

    // Store univariate quadrature in first direction
    forAll(weights1D, wi)
    {
        weights1D[wi] = univariateInverter_().weights()[wi];
        abscissae1D[wi] = univariateInverter_().abscissae()[wi];
    }
}

void Foam::multivariateMomentInversions::CHyQMOMPlus::invert2D
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
    scalar s10 = moments(1)/m00;
    scalar s01 = moments(0, 1)/m00;
    scalar s20 = moments(2)/m00;
    scalar s11 = moments(1, 1)/m00;
    scalar s02 = moments(0, 2)/m00;
    scalar s30 = moments(3)/m00;
    scalar s03 = moments(0, 3)/m00;
    scalar s21 = moments(2, 1)/m00;
    scalar s12 = moments(1, 2)/m00;
    scalar s40 = moments(4)/m00;
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
    scalar c21 = s21;
    scalar c12 = s12;
    scalar c40 = s40;
    scalar c04 = s04;

    // NOTE: Sign changed due to -= operator
    c30 -= (3.0*meanU*s20 - 2.0*pow3(meanU));
    c03 -= (3.0*meanV*s02 - 2.0*pow3(meanV));
    c21 -= (meanV*s20 + 2.0*meanU*s11 - 2.0*sqrMeanU*meanV);
    c12 -= (meanU*s02 + 2.0*meanV*s11 - 2.0*sqrMeanV*meanU);
    c40 -= (4.0*meanU*s30 - 6.0*sqrMeanU*s20 + 3.0*sqr(sqrMeanU));
    c04 -= (4.0*meanV*s03 - 6.0*sqrMeanV*s02 + 3.0*sqr(sqrMeanV));

    realizabilityUnivariateMoments(c20, c30, c40);
    realizabilityUnivariateMoments(c02, c03, c04);

    if (c20 < varMin_)
    {
        univariateMomentSet mDir2({scalar(1), scalar(0), c02, c03, c04}, "R");

        //NOTE: Leave Vf elements null. AP
        univariateInverter_().invert(mDir2);

        scalarList wDir2(univariateInverter_().weights());
        scalarList absDir2(univariateInverter_().abscissae());

        forAll(wDir2, i)
        {
            forAll(wDir2, j)
            {
                if (i == 1)
                {
                    weights2D(i, j) = m00*wDir2[j];
                    abscissae2D(i, j) = vector2D(meanU, absDir2[j] + meanV);
                }
                else
                {
                    weights2D(i, j) = Zero;
                    abscissae2D(i, j) = Zero;
                }
            }
        }

        return;
    }

    // One-dimensional inversion with realizability test
    univariateMomentSet mDir1({scalar(1), scalar(0), c20, c30, c40}, "R");

    // Find univariate quadrature in first direction
    univariateInverter_().invert(mDir1);

    // Store univariate quadrature in first direction
    scalarList wDir1(univariateInverter_().weights());
    scalarList absDir1(univariateInverter_().abscissae());

    if (c02 < varMin_)
    {
        forAll(wDir1, i)
        {
            forAll(wDir1, j)
            {
                if (j == 1)
                {
                    weights2D(i, j) = m00*wDir1[i];
                    abscissae2D(i, j) = vector2D(absDir1[i] + meanU, meanV);
                }
                else
                {
                    weights2D(i, j) = Zero;
                    abscissae2D(i, j) = Zero;
                }
            }
        }

        return;
    }

    // X and y directions are non-degenerate

    // Order of polynomial
    label pOrder = 2;
    scalar sqrtC20 = sqrt(c20);
    scalarList us(absDir1);
    forAll(us, i)
    {
        us[i] /= sqrtC20;
    }

    scalar c11s = c11/sqrtC20;
    scalar c21s = c21/c20;

    scalar q = c30/pow3(sqrtC20);
    scalar eta = c40/sqr(c20);

    // Check for perfect correlation (v = a*u)
    if (sqr(c11s) > c02*(1.0 - SMALL))
    {
        c11s = sign(c11s)*sqrt(c02);
        pOrder = 1;
    }

    scalar r = eta - 1.0 - sqr(q);
    scalar a0 = 0.0;
    scalar a1 = 0.0;

    if (r > minCorrelation_ && pOrder == 2)
    {
        a0 = (c21s - q*c11s)/r;
        a1 = ((eta - 1.0)*c11s - q*c21s)/r;
    }
    else
    {
        a1 = c11s;
    }

    scalarList Vf(3, Zero);
    forAll(Vf, i)
    {
        Vf[i] = a1*us[i] + a0*(sqr(us[i]) - 1.0);
    }

    // Compute conditional variance
    scalar b0 = c02;

    forAll(Vf, vi)
    {
        b0 -= wDir1[vi]*sqr(Vf[vi]);
    }
    scalar b1 = 0.0;

    if (b0 <= 0)
    {
        if (pOrder == 2)
        {
            forAll(Vf, i)
            {
                Vf[i] = c11s*us[i];
            }

            b0 = c02;
            forAll(Vf, vi)
            {
                b0 -= wDir1[vi]*sqr(Vf[vi]);
            }
            pOrder = 1;
        }
    }

    b0 = max(b0, scalar(0));

    if (pOrder == 2)
    {
        b1 = c12/sqrtC20;

        forAll(us, i)
        {
            b1 -= wDir1[i]*sqr(Vf[i])*us[i];
        }
    }

    scalarList mu2(3, b0);

    forAll(mu2, i)
    {
        mu2[i] += b1*us[i];
    }

    label minMu2i = findMin(mu2);
    scalar minMu2 = mu2[minMu2i];

    if (minMu2 < 0)
    {
        b1 = -b0/us[minMu2i];

        forAll(mu2, i)
        {
            mu2[i] = b0 + b1*us[i];
        }

        mu2[minMu2i] = Zero;
    }

    // Check realizability of 3rd and 4th order moments
    q = 0.0;
    eta = 1.0;
    scalar sum1 = 0.0;

    forAll(mu2, i)
    {
        sum1 += wDir1[i]*sqrt(pow3(mu2[i]));
    }

    if (sum1 > SMALL)
    {
        scalar sum03 = c03;

        forAll(Vf, i)
        {
            sum03 -= (wDir1[i]*(pow3(Vf[i]) + 3.0*Vf[i]*mu2[i]));
        }
        q = sum03/sum1;
    }

    scalar sum2 = 0.0;

    forAll(mu2, i)
    {
        sum2 += wDir1[i]*sqr(mu2[i]);
    }

    if (sum1 > SMALL)
    {
        scalar sum04 = c04;

        forAll(Vf, i)
        {
            sum04 -=
                wDir1[i]
                *(
                    pow4(Vf[i])
                  + 6.0*sqr(Vf[i])*mu2[i]
                  + q*4.0*Vf[i]*sqrt(pow3(mu2[i]))
                );
        }

        eta = sum04/sum2;

        if (eta < (sqr(q) + 1.0))
        {
            q = calcQ(q, eta);
            eta = sqr(q) + 1.0;
        }
    }
    scalarList mu3(3, Zero);
    scalarList mu4(3, Zero);

    forAll(mu3, i)
    {
        mu3[i] = q*sqrt(pow3(mu2[i]));
        mu4[i] = eta*sqr(mu2[i]);
    }

    for (label i = 0; i < 3; i++)
    {
        univariateMomentSet mMu({scalar(1), scalar(0), mu2[i], mu3[i], mu4[i]}, "R");
        univariateInverter_().invert(mMu);

        for (label j = 0; j < 3; j++)
        {
            weights2D(i, j) = m00*wDir1[i]*univariateInverter_().weights()[j];

            abscissae2D(i, j) =
                vector2D
                (
                    absDir1[i] + meanU,
                    Vf[i] + univariateInverter_().abscissae()[j] + meanV
                );
        }
    }
}

void Foam::multivariateMomentInversions::CHyQMOMPlus::invert3D
(
    const multivariateMomentSet& moments
)
{
    scalar m000 = moments(0, 0, 0);

    if (m000 < SMALL)
    {
        weights_(1,1,1) = m000;
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
    scalar s210 = moments(2, 1, 0)/m000;
    scalar s201 = moments(2, 0, 1)/m000;
    scalar s120 = moments(1, 2)/m000;
    scalar s111 = moments(1, 1, 1)/m000;
    scalar s102 = moments(1, 0, 2)/m000;
    scalar s030 = moments(0, 3)/m000;
    scalar s021 = moments(0, 2, 1)/m000;
    scalar s012 = moments(0, 1, 2)/m000;
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

    scalar c200 = s200 - sqrMeanU;
    scalar c110 = s110 - meanU*meanV;
    scalar c101 = s101 - meanU*meanW;
    scalar c020 = s020 - sqrMeanV;
    scalar c011 = s011 - meanV*meanW;
    scalar c002 = s002 - sqrMeanW;
    scalar c300 = s300;
    scalar c210 = s210;
    scalar c201 = s201;
    scalar c120 = s120;
    scalar c111 = s111;
    scalar c102 = s102;
    scalar c030 = s030;
    scalar c021 = s021;
    scalar c012 = s012;
    scalar c003 = s003;
    scalar c400 = s400;
    scalar c040 = s040;
    scalar c004 = s004;

    c300 -= (3.0*meanU*s200 - 2.0*pow3(meanU));
    c210 -= (meanV*s200 + 2.0*meanU*s110 - 2.0*sqrMeanU*meanV);
    c201 -= (meanW*s200 + 2.0*meanU*s101 - 2.0*sqrMeanU*meanW);
    c120 -= (meanU*s020 + 2.0*meanV*s110 - 2.0*sqrMeanV*meanU);
    c111 -= (meanU*s011 + meanV*s101 + meanW*s110 - 2.0*meanU*meanV*meanW);
    c102 -= (meanU*s002 + 2.0*meanW*s101 - 2.0*sqrMeanW*meanU);
    c030 -= (3.0*meanV*s020 - 2.0*pow3(meanV));
    c021 -= (meanW*s020 + 2.0*meanV*s011 - 2.0*sqrMeanV*meanW);
    c012 -= (meanV*s002 + 2.0*meanW*s011 - 2.0*sqrMeanW*meanV);
    c003 -= (3.0*meanW*s002 - 2.0*pow3(meanW));

    c400 -= (4.0*meanU*s300 - 6.0*sqrMeanU*s200 + 3.0*sqr(sqrMeanU));
    c040 -= (4.0*meanV*s030 - 6.0*sqrMeanV*s020 + 3.0*sqr(sqrMeanV));
    c004 -= (4.0*meanW*s003 - 6.0*sqrMeanW*s002 + 3.0*sqr(sqrMeanW));

    realizabilityUnivariateMoments(c200, c300, c400);
    realizabilityUnivariateMoments(c020, c030, c040);
    realizabilityUnivariateMoments(c002, c003, c004);

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
                    c021,
                    c012,
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

    // Invert first direction
    univariateMomentSet mDir1({scalar(1), scalar(0), c200, c300, c400}, "R");

    // Find univariate quadrature in first direction
    univariateInverter_().invert(mDir1);

    // Store univariate quadrature in first direction
    scalarList wDir1(univariateInverter_().weights());
    scalarList absDir1(univariateInverter_().abscissae());

    // Y direction is degenerate
    if (c020 < varMin_)
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
                c201,
                c102,
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
            vector2D::zero
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
                c210,
                c120,
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
            vector2D::zero
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
                    velocityAbscissae_(i, j, 1) =
                        vector
                        (
                            abscissaeDir12(i, j).x() + meanU,
                            abscissaeDir12(i, j).y() + meanV,
                            meanW
                        );
                }
            }

            return;
        }
        // All directions are non-degenerate
        label NB = 6;

        for (label i = 0; i < 3; i++)
        {
            wDir12(0, i) /= wDir1[0];
            wDir12(1, i) /= wDir1[1];
            wDir12(2, i) /= wDir1[2];
        }

        // Compute Vf reconstruction
        scalarList Vf(3, Zero);
        for (label i = 0; i < 3; i++)
        {
            Vf[0] += wDir12(0, i)*abscissaeDir12(0, i).y();
            Vf[1] += wDir12(1, i)*abscissaeDir12(1, i).y();
            Vf[2] += wDir12(2, i)*abscissaeDir12(2, i).y();
        }

        mappedList<scalar> absDir2(9, twoDimNodeIndexes, Zero);
        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                absDir2(i, j) = abscissaeDir12(i, j).y() - Vf[i];
            }
        }

        scalar sqrtC200 = sqrt(c200);
        scalar sqrtC020 = sqrt(c020);
        scalar sqrtC002 = sqrt(c002);

        scalarSquareMatrix RAB(3, Zero);
        scalarSquareMatrix Vps(3, Zero);

        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                RAB(i, j) = wDir12(i, j)*wDir1[i];
                Vps(i, j) = absDir2(i, j)/sqrtC020;
            }
        }

        scalarSquareMatrix UABs(3, Zero);
        scalarSquareMatrix VABs(Vps);

        scalarSquareMatrix C00(RAB);
        scalarSquareMatrix C10(RAB);
        scalarSquareMatrix C01(RAB);
        scalarSquareMatrix C11(RAB);
        scalarSquareMatrix C20(RAB);
        scalarSquareMatrix C02(RAB);

        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                UABs(i, j) = absDir1[i]/sqrtC200;
                VABs(i, j) += Vf[i]/sqrtC020;

                C10(i, j) *= UABs(i, j);
                C01(i, j) *= VABs(i, j);
                C11(i, j) *= UABs(i, j)*VABs(i, j);
                C20(i, j) *= sqr(UABs(i, j));
                C02(i, j) *= sqr(VABs(i, j));
            }
        }

        scalarSquareMatrix A(6, 0.0);
        scalarSquareMatrix Vc0(3, 1.0);
        scalarSquareMatrix Vc1(3, 0.0);
        scalarSquareMatrix Vc2(3, 0.0);
        scalarSquareMatrix Vc3(3, 0.0);
        scalarSquareMatrix Vc4(3, 0.0);
        scalarSquareMatrix Vc5(3, 0.0);

        A(0, 0) = 1.0;
        A(0, 4) = 1.0;
        A(1, 1) = 1.0;
        A(4, 0) = 1.0;
        A(5, 0) = 1.0;

        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                Vc1(i, j) = UABs(i, j);
                Vc2(i, j) = Vps(i, j);
                Vc3(i, j) = Vc1(i, j)*Vc2(i, j);
                Vc4(i, j) = sqr(Vc1(i, j));
                Vc5(i, j) = sqr(Vc2(i, j));

                A(0, 5) += C00(i, j)*Vc5(i, j);

                A(1, 4) += C10(i, j)*Vc4(i, j);
                A(1, 5) += C10(i, j)*Vc5(i, j);

                A(2, 1) += C01(i, j)*Vc1(i, j);
                A(2, 2) += C01(i, j)*Vc2(i, j);
                A(2, 3) += C01(i, j)*Vc3(i, j);
                A(2, 4) += C01(i, j)*Vc4(i, j);
                A(2, 5) += C01(i, j)*Vc5(i, j);

                A(3, 0) += C11(i, j);
                A(3, 1) += C11(i, j)*Vc1(i, j);
                A(3, 2) += C11(i, j)*Vc2(i, j);
                A(3, 3) += C11(i, j)*Vc3(i, j);
                A(3, 4) += C11(i, j)*Vc4(i, j);
                A(3, 5) += C11(i, j)*Vc5(i, j);

                A(4, 1) += C20(i, j)*Vc1(i, j);
                A(4, 4) += C20(i, j)*Vc4(i, j);
                A(4, 5) += C20(i, j)*Vc5(i, j);

                A(5, 1) += C02(i, j)*Vc1(i, j);
                A(5, 2) += C02(i, j)*Vc2(i, j);
                A(5, 3) += C02(i, j)*Vc3(i, j);
                A(5, 4) += C02(i, j)*Vc4(i, j);
                A(5, 5) += C02(i, j)*Vc5(i, j);
            }
        }

        scalar c101s = c101/sqrtC200;
        scalar c011s = c011/sqrtC020;
        scalar c111s = c111/sqrtC200/sqrtC020;
        scalar c110s = c110/sqrtC200/sqrtC020;
        scalar c201s = c201/c200;
        scalar c021s = c021/c020;

        if (sqr(c101s) >= c002*(1.0 - SMALL))
        {
            c101s = sign(c101s)*sqrtC002;
            NB = 2;
        }
        else if (sqr(c011s) >= c002*(1.0 - SMALL))
        {
            c011s = sign(c011s)*sqrtC002;
            c101s = c110s*c011s;
            NB = 3;
        }

        scalarRectangularMatrix r(6, 1, 0.0);
        r(1, 0) = c101s;
        r(2, 0) = c011s;
        r(3, 0) = c111s;
        r(4, 0) = c201s;
        r(5, 0) = c021s;
        Foam::SVD svd(A, 1e-3);
        scalarField c(6, Zero);
        scalarRectangularMatrix tmpc(svd.VSinvUt()*r);

        forAll(c, i)
        {
            c[i] = tmpc(i,0);
        }

        scalarSquareMatrix Wf
        (
            c[0]*Vc0 + c[1]*Vc1 + c[2]*Vc2 + c[3]*Vc3 + c[4]*Vc4 + c[5]*Vc5
        );

        scalar sum002 = 0.0;
        scalar sum102 = 0.0;
        scalar sum012 = 0.0;

        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                sum002 += RAB(i, j)*sqr(Wf(i, j));
                sum102 += RAB(i, j)*UABs(i, j)*sqr(Wf(i, j));
                sum012 += RAB(i, j)*VABs(i, j)*sqr(Wf(i, j));
            }
        }

        scalar c102s = c102/sqrtC200;
        scalar c012s = c012/sqrtC020;
        scalar b0 = 1.0;
        forAll(Vf, i)
        {
            b0 -= wDir1[i]*sqr(Vf[i])/c020;
        }
        scalarList d(3, Zero);
        d[0] = c002 - sum002;

        if (NB > 3)
        {
            d[1] = c102s - sum102;
            if (mag(b0) > minCorrelation_)
            {
                d[2] = (c012s - sum012 - d[0]*c110s)/b0;
            }
        }

        scalarSquareMatrix mu2(3, Zero);

        if (d[0] > 0)
        {
            mu2 = d[0]*scalarSquareMatrix(3, 1.0) + d[1]*UABs + d[2]*Vps;
        }

        if (min(mu2) < 0)
        {
            scalarSquareMatrix X(d[1]*UABs + d[2]*Vps);
            scalar x = min(X);
            scalar y = -d[0]/(x - 1e-10);
            mu2 =
                d[0]*scalarSquareMatrix(3, 1.0)
              + y*(d[1]*UABs + d[2]*Vps);
        }

        // Check realizability of 3rd and 4th order moments
        scalar q = 0.0;
        scalar eta = 1.0;
        scalar sum1 = 0.0;
        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                if (mu2(i, j) < 0)
                {
                    mu2(i, j) = 0.0;
                }

                sum1 += RAB(i, j)*sqrt(pow3(mu2(i, j)));
            }
        }

        if (sum1 > SMALL)
        {
            scalar sum3 = c003;
            for (label i = 0; i < 3; i++)
            {
                for (label j = 0; j < 3; j++)
                {
                    sum3 -=
                        RAB(i, j)*(pow3(Wf(i, j)) + 3.0*Wf(i, j)*mu2(i, j));
                }
            }

            q = sum3/sum1;
        }

        scalar sum2 = 0.0;

        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                sum2 += RAB(i, j)*sqr(mu2(i, j));
            }
        }

        if (sum1 > SMALL)
        {
            scalar sum04A = c004;
            scalar sum04B = 0.0;

            for (label i = 0; i < 3; i++)
            {
                for (label j = 0; j < 3; j++)
                {
                    sum04A -=
                        RAB(i, j)
                        *(
                            pow4(Wf(i, j))
                          + 6.0*sqr(Wf(i, j))*mu2(i, j)
                        );

                    sum04B -= 4.0*RAB(i, j)*Wf(i, j)*sqrt(pow3(mu2(i, j)));
                }
            }

            scalar etaA = sum04A/sum2;
            scalar etaB = sum04B/sum2;
            eta = etaA + q*etaB;

            if (eta < (sqr(q) + 1.0))
            {
                q = calcQ(q, eta);
                eta = sqr(q) + 1.0;
            }
        }

        scalarSquareMatrix mu3(3, 0.0);
        scalarSquareMatrix mu4(3, 0.0);

        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                mu3(i, j) = q*sqrt(pow3(mu2(i, j)));
                mu4(i, j) = eta*sqr(mu2(i, j));
            }
        }

        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                univariateMomentSet mMu
                (
                    {scalar(1), scalar(0), mu2(i, j), mu3(i, j), mu4(i, j)},
                    "R"
                );

                univariateInverter_().invert(mMu);

                for (label k = 0; k < 3; k++)
                {
                    weights_(i, j, k) =
                        m000
                       *wDir1[i]
                       *wDir12(i, j)
                       *univariateInverter_().weights()[k];

                    velocityAbscissae_(i, j, k) =
                        vector
                        (
                            absDir1[i] + meanU,
                            Vf[i] + absDir2(i, j) + meanV,
                            Wf(i, j)
                          + univariateInverter_().abscissae()[k]
                          + meanW
                        );
                }
            }
        }
    }
}

bool Foam::multivariateMomentInversions::CHyQMOMPlus::invert
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
