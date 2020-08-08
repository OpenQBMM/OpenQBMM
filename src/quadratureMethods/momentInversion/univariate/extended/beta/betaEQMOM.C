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
    Copyright (C) 2019-2020 Alberto Passalacqua
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

#include "betaEQMOM.H"
#include "scalar.H"
#include "scalarMatrices.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(betaEQMOM, 0);

    addToRunTimeSelectionTable
    (
        extendedMomentInversion,
        betaEQMOM,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::betaEQMOM::betaEQMOM
(
    const dictionary& dict,
    const label nMoments,
    const label nSecondaryNodes
)
:
    extendedMomentInversion(dict, nMoments, nSecondaryNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::betaEQMOM::~betaEQMOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::betaEQMOM::secondaryAbscissa
(
    scalar primaryAbscissa,
    scalar secondaryAbscissa,
    scalar sigma
)
{
    return (1.0 + secondaryAbscissa)/2.0;
}

void Foam::betaEQMOM::momentsStarToMoments
(
    scalar sigma,
    univariateMomentSet& moments,
    const univariateMomentSet& momentsStar
)
{
    label nMom = moments.size();

    if (nMom >= 12)
    {
        FatalErrorInFunction
            << "Moment transformation not implemented."
            << abort(FatalError);
    }

    scalar factor = 1.0 + sigma;

    moments[0] = momentsStar[0];
    moments[1] = momentsStar[1];
    moments[2] = (momentsStar[2] + sigma*momentsStar[1])/factor;

    if (nMom >= 5)
    {
        factor *= 1.0 + 2.0*sigma;
        moments[3] = (momentsStar[3] + sigma*(3.0*momentsStar[2]
                + 2.0*sigma*momentsStar[1])
                )/factor;

        factor *= 1.0 + 3.0*sigma;
        moments[4] = (momentsStar[4] + sigma*(6.0*momentsStar[3]
                + sigma*(11.0*momentsStar[2] + 6.0*momentsStar[1]*sigma))
                )/factor;
    }

    if (nMom >= 7)
    {
        factor *= 1.0 + 4.0*sigma;
        moments[5] = (momentsStar[5] + sigma*(10.0*momentsStar[4]
                + sigma*(35.0*momentsStar[3] + sigma*(50.0*momentsStar[2]
                + 24.0*momentsStar[1]*sigma)))
                )/factor;

        factor *= 1.0 + 5.0*sigma;
        moments[6] = (momentsStar[6] + sigma*(15.0*momentsStar[5]
                + sigma*(85.0*momentsStar[4] + sigma*(225.0*momentsStar[3]
                + sigma*(274.0*momentsStar[2] + 120.0*momentsStar[1]*sigma))))
                )/factor;
    }

    if (nMom >= 9)
    {
        factor *= 1.0 + 6.0*sigma;
        moments[7] = (momentsStar[7] + sigma*(21.0*momentsStar[6]
                + sigma*(175.0*momentsStar[5] + sigma*(735.0*momentsStar[4]
                + sigma*(1624.0*momentsStar[3] + sigma*(1764.0*momentsStar[2]
                + 720.0*momentsStar[1]*sigma)))))
                )/factor;

        factor *= 1.0 + 7.0*sigma;
        moments[8] = (momentsStar[8] + sigma*(28.0*momentsStar[7]
                + sigma*(322.0*momentsStar[6] + sigma*(1960.0*momentsStar[5]
                + sigma*(6769.0*momentsStar[4] + sigma*(13132.0*momentsStar[3]
                + sigma*(13068.0*momentsStar[2]
                + 5040.0*momentsStar[1]*sigma))))))
                )/factor;
    }

    if (nMom >= 11)
    {
        factor *= 1.0 + 8.0*sigma;
        moments[9] = (momentsStar[9] + sigma*(36.0*momentsStar[8]
                + sigma*(546.0*momentsStar[7] + sigma*(4536.0*momentsStar[6]
                + sigma*(22449.0*momentsStar[5] + sigma*(67284.0*momentsStar[4]
                + sigma*(118124.0*momentsStar[3]
                + sigma*(109584.0*momentsStar[2]
                + 40320.0*momentsStar[1]*sigma)))))))
                )/factor;

        factor *= 1.0 + 9.0*sigma;
        moments[10] = (momentsStar[10] + sigma*(45.0*momentsStar[9]
                + sigma*(870.0*momentsStar[8]  + sigma*(9450.0*momentsStar[7]
                + sigma*(63273.0*momentsStar[6]
                + sigma*(269325.0*momentsStar[5]
                + sigma*(723680.0*momentsStar[4]
                + sigma*(1172700.0*momentsStar[3]
                + sigma*(1026576.0*momentsStar[2]
                + 362880.0*momentsStar[1]*sigma))))))))
                )/factor;
    }
}

void Foam::betaEQMOM::momentsToMomentsStar
(
    scalar sigma,
    const univariateMomentSet& moments,
    univariateMomentSet& momentsStar
)
{
    label nMom = moments.size();

    if (nMom >= 12)
    {
        FatalErrorInFunction
            << "The number of moments is too large. The maximum number of"
            << "moments allowed with the beta kernel density function is 11."
            << "Moment transformation not implemented."
            << abort(FatalError);
    }

    momentsStar[0] = moments[0];
    momentsStar[1] = moments[1];
    momentsStar[2] = (1.0 + sigma)*moments[2] - sigma*moments[1];

    if (nMom >= 5)
    {
        momentsStar[3] = (1.0 + 2.0*sigma)*(1.0 + sigma)*moments[3]
                + sigma*(-3.0*(1.0 + sigma)*moments[2]
                + moments[1]*sigma);

        momentsStar[4] = (1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)
                *moments[4] + sigma*(-6.0*(1.0 + 2.0*sigma)*(1.0 + sigma)
                *moments[3] + sigma*(7.0*(1.0 + sigma)*moments[2]
                - moments[1]*sigma));
    }

    if (nMom >= 7)
    {
        momentsStar[5] = (1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[5] + sigma*(-10.0*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[4]
                + sigma*(25.0*(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[3]
                + sigma*(-15.0*(1.0 + sigma)*moments[2] + moments[1]*sigma)));

        momentsStar[6] = (1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[6]
                + sigma*(-15.0*(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[5]
                + sigma*(65.0*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)
                *moments[4] + sigma*(-90.0*(1.0 + 2.0*sigma)*(1.0 + sigma)
                *moments[3] + sigma*(31.0*(1.0 + sigma)*moments[2]
                - moments[1]*sigma))));
    }

    if (nMom >= 9)
    {
        momentsStar[7] = (1.0 + 6.0*sigma)*(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)
                *(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[7]
                + sigma*(-21.0*(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)
                *(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[6]
                + sigma*(140.0*(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[5]
                + sigma*(-350.0*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[4] + sigma*(301.0*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[3] + sigma*(-63.0*(1.0 + sigma)
                *moments[2] + moments[1]*sigma)))));

        momentsStar[8] = (1.0 + 7.0*sigma)*(1.0 + 6.0*sigma)*(1.0 + 5.0*sigma)
                *(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[8] + sigma*(-28.0*(1.0 + 6.0*sigma)
                *(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[7]
                + sigma*(266.0*(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)
                *(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[6]
                + sigma*(-1050.0*(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[5]
                + sigma*(1701.0*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[4] + sigma*(-966.0*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[3] + sigma*(127.0*(1.0 + sigma)
                *moments[2] - moments[1]*sigma))))));
    }

    if (nMom >= 11)
    {
        momentsStar[9] = (1.0 + 8.0*sigma)*(1.0 + 7.0*sigma)*(1.0 + 6.0*sigma)
                *(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[9]
                + sigma*(-36.0*(1.0 + 7.0*sigma)*(1.0 + 6.0*sigma)
                *(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[8]
                + sigma*(462.0*(1.0 + 6.0*sigma)*(1.0 + 5.0*sigma)
                *(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[7] + sigma*(-2646.0*(1.0 + 5.0*sigma)
                *(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[6] + sigma*(6951.0*(1.0 + 4.0*sigma)
                *(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[5]
                + sigma*(-7770.0*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[4] + sigma*(3025.0*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[3] + sigma*(-255.0*(1.0 + sigma)
                *moments[2] + moments[1]*sigma)))))));

        momentsStar[10] = (1.0 + 9.0*sigma)*(1.0 + 8.0*sigma)*(1.0 + 7.0*sigma)
                *(1.0 + 6.0*sigma)*(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)
                *(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[10]
                + sigma*(-45.0*(1.0 + 8.0*sigma)*(1.0 + 7.0*sigma)
                *(1.0 + 6.0*sigma)*(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)
                *(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[9]
                + sigma*(750.0*(1.0 + 7.0*sigma)*(1.0 + 6.0*sigma)
                *(1.0 + 5.0*sigma)*(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)
                *(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[8]
                + sigma*(-5880.0*(1.0 + 6.0*sigma)*(1.0 + 5.0*sigma)
                *(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[7] + sigma*(22827.0*(1.0 + 5.0*sigma)
                *(1.0 + 4.0*sigma)*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[6] + sigma*(-42525.0*(1.0 + 4.0*sigma)
                *(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)*(1.0 + sigma)*moments[5]
                + sigma*(34105.0*(1.0 + 3.0*sigma)*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[4] + sigma*(-9330.0*(1.0 + 2.0*sigma)
                *(1.0 + sigma)*moments[3] + sigma*(511.0*(1.0 + sigma)
                *moments[2] - moments[1]*sigma))))))));
    }
}

Foam::scalar Foam::betaEQMOM::m2N
(
    scalar sigma,
    const univariateMomentSet& momentsStar
)
{
    univariateMomentSet mStar(momentsStar);
    label nMomentsStar = mStar.size();

    if (mStar.nRealizableMoments() >= nMomentsStar - 1)
    {
        univariateMomentSet m(nMomentsStar, "01");
        momentsStarToMoments(sigma, m, mStar);

        return m.last();
    }

    return GREAT;
}

void Foam::betaEQMOM::recurrenceRelation
(
    scalarList& a,
    scalarList& b,
    scalar primaryAbscissa,
    scalar sigma
)
{
    scalar alpha = (1.0 - primaryAbscissa)/sigma - 1.0;
    scalar beta = primaryAbscissa/sigma - 1.0;
    scalar nab = 0.0;

    a[0] = (beta - alpha)/(alpha + beta + 2.0);

    b[0] = exp((alpha + beta + 1.0)*log(2.0) + lgamma(alpha + 1.0)
            + lgamma(beta + 1.0) - lgamma(alpha + beta + 2.0));

    a[1] = (sqr(beta) - sqr(alpha))/((alpha + beta + 2.0)*(alpha + beta + 4.0));

    b[1] = (4.0*(alpha + 1.0)*(beta + 1.0))
            /(sqr(alpha + beta + 2.0)*(alpha + beta + 3.0));

    for (label i = 2; i < a.size(); i++)
    {
        nab = 2.0*scalar(i) + alpha + beta;

        a[i] = (sqr(beta) - sqr(alpha))/(nab*(nab + 2.0));

        b[i] = (4.0*(scalar(i) + alpha)*(scalar(i) + beta)*scalar(i)
                *(scalar(i)+ alpha + beta))/(sqr(nab)*(nab + 1.0)*(nab - 1.0));
    }
}

Foam::scalar Foam::betaEQMOM::sigmaMax(univariateMomentSet& moments)
{
    scalar p2 = (moments[0]*moments[2] - sqr(moments[1]))
                /(moments[0]*moments[1] - sqr(moments[1]));

    //Info << "p2 = " << p2 << endl;

    return p2/(1.0 - p2);
}

Foam::tmp<Foam::scalarField> Foam::betaEQMOM::f(const scalarField& x) const
{
    tmp<scalarField> tmpY
    (
        new scalarField(x.size(), Zero)
    );
    scalarField& y = tmpY.ref();

    for (label pNodei = 0; pNodei < nPrimaryNodes_; pNodei++)
    {
        scalar pAbscissa = primaryAbscissae_[pNodei];
        scalar lambda = pAbscissa/sigma_;
        scalar mu = (1.0 - pAbscissa)/sigma_;

        y +=
            pow(x, lambda - 1)*pow(1.0 - x, mu - 1)
           /(tgamma(lambda)*tgamma(mu)/tgamma(lambda + mu))
           *this->primaryWeights_[pNodei];
    }

    return tmpY;
}

// ************************************************************************* //
