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

                // Remove small negative values in abscissae
                scalar bAbscissa1 = max(pAbscissa1[celli], 0.0);
                scalar bAbscissa2 = max(pAbscissa2[celli], 0.0);

                aSource +=
                    pWeight1[celli]*
                    (
                        pWeight2[celli]*
                        (
                            0.5*pow // Birth
                            (
                                pow3(bAbscissa1) + pow3(bAbscissa2),
                                momentOrder/3.0
                            )
                            - pow(bAbscissa1, momentOrder)
                        )*aggregationKernel_->Ka
                            (
                                bAbscissa1, bAbscissa2, celli
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

                    // Remove small negative values in abscissae
                    scalar bAbscissa1 = max(sAbscissa1[celli], 0.0);
                    scalar bAbscissa2 = max(sAbscissa2[celli], 0.0);

                    aSource +=
                        pWeight1[celli]*sWeight1[celli]*
                        (
                            pWeight2[celli]*sWeight2[celli]*
                            (
                                0.5*pow // Birth
                                (
                                    pow3(bAbscissa1) + pow3(bAbscissa2),
                                    momentOrder/3.0
                                )
                              - pow(bAbscissa1, momentOrder)
                            )*aggregationKernel_->Ka
                                (
                                    bAbscissa1, bAbscissa2, celli
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

            scalar bAbscissa = max(node.primaryAbscissa()[celli], 0.0);

            bSource += node.primaryWeight()[celli]
                    *breakupKernel_->Kb(bAbscissa, celli)
                    *(
                        daughterDistribution_->mD(momentOrder, bAbscissa)//Birth
                      - pow(bAbscissa, momentOrder)   //Death
                    );
        }

        return bSource;
    }

    forAll(nodes, pNodeI)
    {
        const volScalarNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            scalar bAbscissa
                = max(node.secondaryAbscissae()[sNodei][celli], 0.0);

            bSource += node.primaryWeight()[celli]
                *node.secondaryWeights()[sNodei][celli]
                *breakupKernel_->Kb(bAbscissa, celli)
                *(
                    daughterDistribution_->mD(momentOrder, bAbscissa)   //Birth
                  - pow(bAbscissa, momentOrder)                         //Death
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

            scalar bAbscissa = max(node.primaryAbscissa()[celli], 0.0);

            gSource += node.primaryWeight()[celli]
                    *growthModel_->Kg(node.primaryAbscissa()[celli])
                    *momentOrder*pow(bAbscissa, momentOrder - 1);
        }

        return gSource;
    }

    forAll(nodes, pNodeI)
    {
        const volScalarNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            scalar bAbscissa
                = max(node.secondaryAbscissae()[sNodei][celli], 0.0);

            gSource += node.primaryWeight()[celli]
                *node.secondaryWeights()[sNodei][celli]
                *growthModel_->Kg(bAbscissa)
                *momentOrder*pow
                    (
                        bAbscissa,
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

    return impSource;
}

Foam::scalar 
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::cellMomentSource
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

Foam::scalar 
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::realizableCo() const
{
    return univariatePDFTransportModel::realizableCo();
}

Foam::scalar 
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::CoNum() const
{
    return 0.0;
}

void 
Foam::PDFTransportModels::populationBalanceModels::univariatePopulationBalance
::solve()
{
    univariatePDFTransportModel::solve();
}

// ************************************************************************* //
