/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

#include "univariatePopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(univariatePopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        univariatePopulationBalance,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::univariatePopulationBalance
(
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    univariatePDFTransportModel(dict, U.mesh(), U, "RPlus"),
    populationBalanceModel(dict, U, phi),
    aggregation_(dict.lookup("aggregation")),
    breakup_(dict.lookup("breakup")),
    growth_(dict.lookup("growth")),
    aggregationKernel_
    (
        Foam::populationBalanceSubModels::aggregationKernel::New
        (
            dict.subDict("aggregationKernel")
        )
    ),
    breakupKernel_
    (
        Foam::populationBalanceSubModels::breakupKernel::New
        (
            dict.subDict("breakupKernel")
        )
    ),
    daughterDistribution_
    (
        Foam::populationBalanceSubModels::daughterDistribution::New
        (
            dict.subDict("daughterDistribution")
        )
    ),
    growthModel_
    (
        Foam::populationBalanceSubModels::growthModel::New
        (
            dict.subDict("growthModel")
        )
    ),
    diffusionModel_
    (
        Foam::populationBalanceSubModels::diffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::~univariatePopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::aggregationSource
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> aSource
    (
        new volScalarField
        (
            IOobject
            (
                "aSource",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    if (!aggregation_)
    {
        aSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return aSource;
    }

    label order = moment.order();

    volScalarField& aggregationSource = aSource.ref();

    forAll(quadrature_.nodes(), pNode1I)
    {
        const extendedVolScalarNode& node1 = quadrature_.nodes()[pNode1I];

        const volScalarField& pWeight1 = node1.primaryWeight();

        forAll(node1.secondaryWeights(), sNode1I)
        {

            const volScalarField& sWeight1 = node1.secondaryWeights()[sNode1I];

            const volScalarField& sAbscissa1
                = node1.secondaryAbscissae()[sNode1I];

            forAll(quadrature_.nodes(), pNode2I)
            {
                const extendedVolScalarNode& node2 = quadrature_.nodes()[pNode2I];

                const volScalarField& pWeight2 = node2.primaryWeight();

                forAll(node2.secondaryWeights(), sNode2I)
                {
                    const volScalarField& sWeight2
                        = node2.secondaryWeights()[sNode2I];

                    const volScalarField& sAbscissa2
                        = node2.secondaryAbscissae()[sNode2I];

                    tmp<volScalarField> aggInnerSum =
                        pWeight1*sWeight1*
                        (
                            pWeight2*sWeight2*
                            (
                                0.5*pow // Birth
                                (
                                    pow3(sAbscissa1) + pow3(sAbscissa2),
                                    order/3.0
                                )
                              - pow(sAbscissa1, order)
                            )*aggregationKernel_->Ka(sAbscissa1, sAbscissa2)
                        );

                    aggregationSource.dimensions().reset
                    (
                        aggInnerSum().dimensions()
                    );

                    aggregationSource == aggregationSource + aggInnerSum();
                }
            }
        }
    }

    return aSource;
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::breakupSource
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> bSource
    (
        new volScalarField
        (
            IOobject
            (
                "bSource",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    if (!breakup_)
    {
        bSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return bSource;
    }

    label order = moment.order();

    volScalarField& breakupSource = bSource.ref();

    forAll(quadrature_.nodes(), pNodeI)
    {
        const extendedVolScalarNode& node = quadrature_.nodes()[pNodeI];

        forAll(node.secondaryWeights(), sNodeI)
        {
            tmp<volScalarField> bSrc = node.primaryWeight()
                *node.secondaryWeights()[sNodeI]
                *breakupKernel_->Kb(node.secondaryAbscissae()[sNodeI])
                *(
                    daughterDistribution_->mD                      //Birth
                    (
                        order,
                        node.secondaryAbscissae()[sNodeI]
                    )
                  - pow(node.secondaryAbscissae()[sNodeI], order)   //Death
                 );

            breakupSource.dimensions().reset(bSrc().dimensions());
            breakupSource == breakupSource + bSrc;
        }
    }

    return bSource;
}

Foam::tmp<fvScalarMatrix> Foam::PDFTransportModels::populationBalanceModels
::univariatePopulationBalance::momentDiffusion
(
    const volUnivariateMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}

Foam::tmp<Foam::volScalarField>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::phaseSpaceConvection
(
    const volUnivariateMoment& moment
)
{
    tmp<volScalarField> gSource
    (
        new volScalarField
        (
            IOobject
            (
                "gSource",
                U_.mesh().time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            U_.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    if (!growth_)
    {
        gSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return gSource;
    }

    label order = moment.order();

    if (order < 1)
    {
        gSource.ref().dimensions().reset(moment.dimensions()/dimTime);

        return gSource;
    }

    volScalarField& growthSource = gSource.ref();

    forAll(quadrature_.nodes(), pNodeI)
    {
        const extendedVolScalarNode& node = quadrature_.nodes()[pNodeI];

        forAll(node.secondaryWeights(), sNodeI)
        {
            tmp<volScalarField> gSrc = node.primaryWeight()
                *node.secondaryWeights()[sNodeI]
                *growthModel_->Kg(node.secondaryAbscissae()[sNodeI])
                *order*pow(node.secondaryAbscissae()[sNodeI],order-1);

            growthSource.dimensions().reset(gSrc().dimensions());
            growthSource == growthSource + gSrc;
        }
    }

    return gSource;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::momentSource
(
    const volUnivariateMoment& moment
)
{
    tmp<fvScalarMatrix> mSource
    (
        new fvScalarMatrix
        (
            moment,
            moment.dimensions()*dimVol/dimTime
        )
    );

    mSource.ref() += aggregationSource(moment) + breakupSource(moment);

    return mSource;
}

void Foam::PDFTransportModels::populationBalanceModels
::univariatePopulationBalance::solve
()
{
    univariatePDFTransportModel::solve();
}

// ************************************************************************* //
