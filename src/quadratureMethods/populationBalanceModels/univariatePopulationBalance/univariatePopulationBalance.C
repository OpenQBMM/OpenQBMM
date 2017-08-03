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
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    univariatePDFTransportModel(name, dict, phi.mesh(), phi, "RPlus"),
    populationBalanceModel(name, dict, phi),
    name_(name),
    aggregation_(dict.lookup("aggregation")),
    breakup_(dict.lookup("breakup")),
    growth_(dict.lookup("growth")),
    aggregationKernel_
    (
        Foam::populationBalanceSubModels::aggregationKernel::New
        (
            dict.subDict("aggregationKernel"),
            phi_.mesh()
        )
    ),
    breakupKernel_
    (
        Foam::populationBalanceSubModels::breakupKernel::New
        (
            dict.subDict("breakupKernel"),
            phi_.mesh()
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
    ),
    nucleationModel_
    (
        Foam::populationBalanceSubModels::nucleationModel::New
        (
            dict.subDict("nucleationModel"),
            phi_.mesh()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::~univariatePopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::aggregationSource
(
    const label& momentOrder,
    const label& celli
)
{
    scalar aSource = 0.0;

    if (!aggregation_)
    {
        return aSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    if (!nodes[0].extended())   // Non-extended quadrature case
    {
        forAll(nodes, pNode1i)
        {
            const volScalarNode& node1 = nodes[pNode1i];
            const volScalarField& pWeight1 = node1.primaryWeight();
            const volScalarField& pAbscissa1 = node1.primaryAbscissa();

            forAll(nodes, pNode2i)
            {
                const volScalarNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();
                const volScalarField& pAbscissa2 = node2.primaryAbscissa();

                aSource +=
                    pWeight1[celli]*
                    (
                        pWeight2[celli]*
                        (
                            0.5*pow // Birth
                            (
                                pow3(pAbscissa1[celli])
                              + pow3(pAbscissa2[celli]),
                                momentOrder/3.0
                            )
                            - pow(pAbscissa1[celli], momentOrder)
                        )*aggregationKernel_->Ka
                            (
                                pAbscissa1[celli], pAbscissa2[celli], celli
                            )
                    );
            }
        }

        return aSource;
    }

    forAll(nodes, pNode1i)      // Extended quadrature case
    {
        const volScalarNode& node1 = nodes[pNode1i];
        const volScalarField& pWeight1 = node1.primaryWeight();

        forAll(node1.secondaryWeights(), sNode1i)
        {
            const volScalarField& sWeight1 = node1.secondaryWeights()[sNode1i];

            const volScalarField& sAbscissa1
                = node1.secondaryAbscissae()[sNode1i];

            forAll(nodes, pNode2i)
            {
                const volScalarNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();

                forAll(node2.secondaryWeights(), sNode2i)
                {
                    const volScalarField& sWeight2
                        = node2.secondaryWeights()[sNode2i];

                    const volScalarField& sAbscissa2
                        = node2.secondaryAbscissae()[sNode2i];

                    aSource +=
                        pWeight1[celli]*sWeight1[celli]*
                        (
                            pWeight2[celli]*sWeight2[celli]*
                            (
                                0.5*pow // Birth
                                (
                                    pow3(sAbscissa1[celli])
                                  + pow3(sAbscissa2[celli]),
                                    momentOrder/3.0
                                )
                              - pow(sAbscissa1[celli], momentOrder)
                            )*aggregationKernel_->Ka
                                (
                                    sAbscissa1[celli], sAbscissa2[celli], celli
                                )
                        );
                }
            }
        }
    }

    return aSource;
}

Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::breakupSource
(
    const label& momentOrder,
    const label& celli
)
{
    scalar bSource = 0.0;

    if (!breakup_)
    {
        return bSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volScalarNode& node = nodes[pNodeI];

            bSource += node.primaryWeight()[celli]
                    *breakupKernel_->Kb(node.primaryAbscissa()[celli], celli)
                    *(
                        daughterDistribution_->mD                      //Birth
                        (
                            momentOrder,
                            node.primaryAbscissa()[celli]
                        )
                    - pow(node.primaryAbscissa()[celli], momentOrder)   //Death
                    );
        }

        return bSource;
    }

    forAll(nodes, pNodeI)
    {
        const volScalarNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            bSource += node.primaryWeight()[celli]
                *node.secondaryWeights()[sNodei][celli]
                *breakupKernel_->Kb
                    (
                        node.secondaryAbscissae()[sNodei][celli], celli
                    )
                *(
                    daughterDistribution_->mD                      //Birth
                    (
                        momentOrder,
                        node.secondaryAbscissae()[sNodei][celli]
                    )                                               //Death
                  - pow(node.secondaryAbscissae()[sNodei][celli], momentOrder)
                 );
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

Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::phaseSpaceConvection
(
    const label& momentOrder,
    const label& celli
)
{
    scalar gSource = 0.0;

    if (!growth_ || momentOrder < 1)
    {
        return gSource;
    }

    const PtrList<volScalarNode>& nodes = quadrature_.nodes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volScalarNode& node = nodes[pNodeI];

            gSource += node.primaryWeight()[celli]
                    *growthModel_->Kg(node.primaryAbscissa()[celli])
                    *momentOrder*pow
                        (
                            node.primaryAbscissa()[celli], momentOrder - 1
                        );
        }

        return gSource;
    }

    forAll(nodes, pNodeI)
    {
        const volScalarNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            gSource += node.primaryWeight()[celli]
                *node.secondaryWeights()[sNodei][celli]
                *growthModel_->Kg(node.secondaryAbscissae()[sNodei][celli])
                *momentOrder*pow
                    (
                        node.secondaryAbscissae()[sNodei][celli],
                        momentOrder - 1
                    );
        }
    }

    return gSource;
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::implicitMomentSource
(
    const volUnivariateMoment& moment
)
{
    tmp<fvScalarMatrix> impSource
    (
        new fvScalarMatrix
        (
            moment,
            moment.dimensions()*dimVol/dimTime
        )
    );

//     impSource.ref() +=
//         aggregationSource(moment) + breakupSource(moment)
//         + nucleationModel_->nucleationSource(moment);

    return impSource;
}

Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::univariatePopulationBalance::cellMomentSource
(
    label& momentOrder,
    label& celli
)
{
    return aggregationSource(momentOrder, celli)
            + breakupSource(momentOrder, celli)
            + nucleationModel_->nucleationSource(momentOrder, celli)
            + phaseSpaceConvection(momentOrder, celli);
}

Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::univariatePopulationBalance::realizableCo
()
{
    return univariatePDFTransportModel::realizableCo();
}

void Foam::PDFTransportModels::populationBalanceModels
::univariatePopulationBalance::solve
()
{
    univariatePDFTransportModel::solve();
}

// ************************************************************************* //
