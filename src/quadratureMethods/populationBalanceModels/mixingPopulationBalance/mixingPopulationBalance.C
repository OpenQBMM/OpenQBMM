/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
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

#include "mixingPopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(mixingPopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        populationBalanceModel,
        mixingPopulationBalance,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::mixingPopulationBalance
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    populationBalanceModel(name, dict, phi),
    ATol_(readScalar(dict.subDict("odeCoeffs").lookup("ATol"))),
    RTol_(readScalar(dict.subDict("odeCoeffs").lookup("RTol"))),
    fac_(readScalar(dict.subDict("odeCoeffs").lookup("fac"))),
    facMin_(readScalar(dict.subDict("odeCoeffs").lookup("facMin"))),
    facMax_(readScalar(dict.subDict("odeCoeffs").lookup("facMax"))),
    minLocalDt_(readScalar(dict.subDict("odeCoeffs").lookup("minLocalDt"))),
    name_(name),
    mixingModel_
    (
        new Foam::PDFTransportModels::mixingModels::turbulentMixing
        (
            "mixing",
            dict.subDict("mixing"),
            phi
        )
    ),
    minMixtureFractionVariance_
    (
        dict.lookupOrDefault("minMixtureFractionVariance", 1.0e-4)
    ),
    minEnvironmentWeight_
    (
        dict.lookupOrDefault("minEnvironmentWeight", 1.0e-6)
    ),
    p1_
    (
        mixingModel_().quadrature().nodes()[0].primaryWeight()
    ),
    xi1_
    (
        mixingModel_().quadrature().nodes()[0].primaryAbscissa()
    ),
    p2_
    (
        mixingModel_().quadrature().nodes()[1].primaryWeight()
    ),
    xi2_
    (
        mixingModel_().quadrature().nodes()[1].primaryAbscissa()
    ),
    meanXi_
    (
        mixingModel_().quadrature().moments()[1]
    ),
    meanMomentsQuadrature_
    (
        name + "MeanMoments",
        phi_.mesh(),
        "RPlus"
    ),
    meanMomentsVarianceQuadrature_
    (
        name + "MeanMomentsVariance",
        phi_.mesh(),
        "RPlus"
    ),
    meanMomentsAdvection_
    (
        univariateMomentAdvection::New
        (
            meanMomentsQuadrature_.subDict("momentAdvection"),
            meanMomentsQuadrature_,
            phi_,
            "RPlus"
        )
    ),
    meanMomentsVarianceAdvection_
    (
        univariateMomentAdvection::New
        (
            meanMomentsVarianceQuadrature_.subDict("momentAdvection"),
            meanMomentsVarianceQuadrature_,
            phi_,
            "RPlus"
        )
    ),
    meanMoments_
    (
        meanMomentsQuadrature_.moments()
    ),
    meanMomentsVariance_
    (
        meanMomentsVarianceQuadrature_.moments()
    ),
    envOneQuadrature_
    (
        IOobject::groupName
        (
            "quadratureProperties",
            meanMomentsQuadrature_.name()
        ),
        "envOneQuadrature",
        meanMomentsQuadrature_.moments()
    ),
    envTwoQuadrature_
    (
        IOobject::groupName
        (
            "quadratureProperties",
            meanMomentsQuadrature_.name()
        ),
        "envTwoQuadrature",
        meanMomentsQuadrature_.moments()
    ),
    mEnvOne_(envOneQuadrature_.moments()),
    mEnvTwo_(envTwoQuadrature_.moments()),
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
    ),
    envMixingModel_
    (
        Foam::populationBalanceSubModels::environmentMixingModel::New
        (
            dict.subDict("environmentMixingModel"),
            phi_.mesh()
        )
    )

{
    if (mixingModel_().quadrature().nodes().size() != 2)
    {
         FatalErrorInFunction
            << "The mixingPbe model can only be used with two environments."
            << endl << "The mixing model must use two quadrature nodes."
            << abort(FatalError);
    }

    // Compute moments in the environments
    calcEnvironmentMoments();

    // Update quadrature in the environments
    envOneQuadrature_.updateQuadrature();
    envTwoQuadrature_.updateQuadrature();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::~mixingPopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::explicitMomentSource
(
    univariateQuadratureApproximation& quadrature,
    const label environment
)
{
    volUnivariateMomentFieldSet& moments(quadrature.moments());
    const mappedPtrList<volScalarNode>& nodes(quadrature.nodes());
    label nMoments = quadrature.nMoments();
    scalar globalDt = moments[0].mesh().time().deltaT().value();

    Info << "Solving source terms in realizable ODE solver." << endl;

    forAll(moments[0], celli)
    {
        // Storing old moments to recover from failed step
        scalarList oldMoments(nMoments, 0.0);

        forAll(oldMoments, mi)
        {
            oldMoments[mi] = moments[mi][celli];
        }

        //- Local time
        scalar localT = 0.0;

        // Initialize the local step
        scalar localDt = globalDt/100.0;

        // Initialize RK parameters
        scalarList k1(nMoments, 0.0);
        scalarList k2(nMoments, 0.0);
        scalarList k3(nMoments, 0.0);

        // Flag to indicate if the time step is complete
        bool timeComplete = false;

        // Check realizability of intermediate moment sets
        bool realizableUpdate1 = false;
        bool realizableUpdate2 = false;
        bool realizableUpdate3 = false;

        scalarList momentsSecondStep(nMoments, 0.0);

        bool nullSource = false;

        while (!timeComplete)
        {
            do
            {
                // First intermediate update
                forAll(oldMoments, mi)
                {
                    k1[mi]
                        = localDt*cellMomentSource
                          (
                              mi, celli, nodes, environment
                          );

                    moments[mi][celli] = oldMoments[mi] + k1[mi];

                    nullSource = (mag(k1[mi]) < SMALL) || nullSource;
                }

                if (nullSource)
                {
                    return;
                }

                realizableUpdate1 =
                        quadrature.updateLocalQuadrature(celli, false);

                quadrature.updateLocalMoments(celli);

                // Second moment update
                forAll(oldMoments, mi)
                {
                    k2[mi]
                        = localDt*cellMomentSource
                          (
                              mi, celli, nodes, environment
                          );

                    moments[mi][celli] = oldMoments[mi] + (k1[mi] + k2[mi])/4.0;
                    momentsSecondStep[mi] = moments[mi][celli];
                }

                realizableUpdate2 =
                        quadrature.updateLocalQuadrature(celli, false);

                quadrature.updateLocalMoments(celli);

                // Third moment update
                forAll(oldMoments, mi)
                {
                    k3[mi]
                        = localDt*cellMomentSource
                          (
                              mi, celli, nodes, environment
                          );

                    moments[mi][celli] =
                        oldMoments[mi] + (k1[mi] + k2[mi] + 4.0*k3[mi])/6.0;
                }

                realizableUpdate3 =
                        quadrature.updateLocalQuadrature(celli, false);

                quadrature.updateLocalMoments(celli);

                if
                (
                    !realizableUpdate1
                 || !realizableUpdate2
                 || !realizableUpdate3
                )
                {
                    localDt /= 2.0;

                    if (localDt < minLocalDt_)
                    {
                        FatalErrorInFunction
                            << "Reached minimum local step in realizable ODE"
                            << nl
                            << "    solver. Cannot ensure realizability." << nl
                            << "    Local time step = " << localDt << nl
                            << "    Min local time step = " << minLocalDt_ << nl
                            << "    Last valid moments in cell: "
                            << oldMoments << nl
                            << abort(FatalError);
                    }

                    forAll(oldMoments, mi)
                    {
                        moments[mi][celli] = oldMoments[mi];
                    }
                }
            }
            while
            (
                !realizableUpdate1
             || !realizableUpdate2
             || !realizableUpdate3
            );

            scalar error = 0.0;

            for (label mi = 0; mi < nMoments; mi++)
            {
                scalar scalei =
                        ATol_
                    + max
                        (
                            mag(momentsSecondStep[mi]), mag(oldMoments[mi])
                        )*RTol_;

                error +=
                        sqr
                        (
                            (momentsSecondStep[mi] - moments[mi][celli])/scalei
                        );
            }

            error = sqrt(error/nMoments);

            if (error < SMALL)
            {
                return;
            }

            if (error < 1)
            {
                localDt *= min(facMax_, max(facMin_, fac_/pow(error, 1.0/3.0)));

                scalar maxLocalDt = max(globalDt - localT, 0.0);
                localDt = min(maxLocalDt, localDt);

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }

                if (localDt == 0.0)
                {
                    timeComplete = true;
                    localT = 0.0;
                    break;
                }

                localT += localDt;
            }
            else
            {
                localDt *= min(1.0, max(facMin_, fac_/pow(error, 1.0/3.0)));

                forAll(oldMoments, mi)
                {
                    moments[mi][celli] = oldMoments[mi];
                }
            }
        }
    }
}


void Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::calcEnvironmentMoments()
{
    const volUnivariateMomentFieldSet& mXi
    (
        mixingModel_().quadrature().moments()
    );

    const volUnivariateMoment& xiMean_ = mXi[1];
    const volUnivariateMoment& xiMTwo_ = mXi[2];

    // Compute variance of the mixture fraction
    const volScalarField xiVariance(xiMTwo_ - sqr(xiMean_));

    // Difference between abscissae of the mixture fraction
    volScalarField xiDiff(xi1_ - xi2_);

    // Compute moments in the environment
    forAll(xiDiff, celli)
    {
        // Null or very small variance of the mixture fraction
        // Moments in the two environments are identical
        if (xiDiff[celli] > minMixtureFractionVariance_)
        {
            forAll(mEnvOne_, mi)
            {
                if (p1_[celli] > minEnvironmentWeight_)
                {
                    mEnvOne_[mi][celli]
                        =
                        (
                            meanMomentsVariance_[mi][celli]
                           - meanMoments_[mi][celli]*xi2_[celli]
                        )/(p1_[celli]*xiDiff[celli]);
                }
                else
                {
                    mEnvOne_[mi][celli] = 0;
                }

                if (p2_[celli] > minEnvironmentWeight_)
                {
                    mEnvTwo_[mi][celli]
                        =
                        (
                            meanMoments_[mi][celli]*xi1_[celli]
                          - meanMomentsVariance_[mi][celli]
                        )/(p2_[celli]*xiDiff[celli]);
                }
                else
                {
                    mEnvTwo_[mi][celli] = 0;
                }
            }
        }
        else
        {
            forAll(mEnvOne_, mi)
            {
                mEnvOne_[mi][celli] = meanMoments_[mi][celli];
                mEnvTwo_[mi][celli] = meanMoments_[mi][celli];
            }
        }
    }

    forAll(mEnvOne_, mi)
    {
        mEnvOne_[mi].correctBoundaryConditions();
        mEnvTwo_[mi].correctBoundaryConditions();
    }
}

void Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::calcMixedMoments()
{
    forAll(meanMoments_, mi)
    {
        meanMoments_[mi] == p1_*mEnvOne_[mi] + p2_*mEnvTwo_[mi];

        meanMomentsVariance_[mi]
            == p1_*xi1_*mEnvOne_[mi] + p2_*xi2_*mEnvTwo_[mi];
    }
}

Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::aggregationSource
(
    const label momentOrder,
    const label celli,
    const mappedPtrList<volScalarNode>& nodes,
    const label environment
)
{
    scalar aSource = 0.0;

    if (!aggregation_)
    {
        return aSource;
    }

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
                                pAbscissa1[celli],
                                pAbscissa2[celli],
                                celli,
                                environment
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
                                    sAbscissa1[celli],
                                    sAbscissa2[celli],
                                    celli,
                                    environment
                                )
                        );
                }
            }
        }
    }

    return aSource;
}

Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::breakupSource
(
    const label momentOrder,
    const label celli,
    const mappedPtrList<volScalarNode>& nodes,
    const label environment
)
{
    scalar bSource = 0.0;

    if (!breakup_)
    {
        return bSource;
    }

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
::mixingPopulationBalance::momentDiffusion
(
    const volUnivariateMoment& moment
)
{
    return diffusionModel_->momentDiff(moment);
}

Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::phaseSpaceConvection
(
    const label momentOrder,
    const label celli,
    const mappedPtrList<volScalarNode>& nodes,
    const label environment
)
{
    scalar gSource = 0.0;

    if (!growth_ || momentOrder < 1)
    {
        return gSource;
    }

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

Foam::scalar 
Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::cellMomentSource
(
    label momentOrder,
    label celli,
    const mappedPtrList<volScalarNode>& nodes,
    const label environment
)
{
    return aggregationSource(momentOrder, celli, nodes, environment)
            + breakupSource(momentOrder, celli, nodes, environment)
            + nucleationModel_->nucleationSource(momentOrder, celli)
            + phaseSpaceConvection(momentOrder, celli, nodes, environment);
}

Foam::scalar 
Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::realizableCo() const
{
    return
        min
        (
            meanMomentsAdvection_().realizableCo(),
            meanMomentsVarianceAdvection_().realizableCo()
        );
}

Foam::scalar 
Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::CoNum() const
{
    return 0.0;
}

void Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::solve()
{
    // Solve mixingModel
    Info << "Solving mixing model" << endl;

    mixingModel_().solve();

    Info << "Solving population balance" << endl;

    // Update PBE advection
    meanMomentsAdvection_().update();
    meanMomentsVarianceAdvection_().update();

    // List of moment transport equations for the mean moments
    PtrList<fvScalarMatrix> meanMomentEqns
    (
        meanMomentsQuadrature_.nMoments()
    );

    // List of moment transport equations for the variance of the mean moments
    PtrList<fvScalarMatrix> meanMomentVarianceEqns
    (
        meanMomentsVarianceQuadrature_.nMoments()
    );

    // Solve moment transport equations
    forAll(meanMomentsQuadrature_.moments(), momenti)
    {
        volUnivariateMoment& meanM
        (
            meanMomentsQuadrature_.moments()[momenti]
        );

        meanMomentEqns.set
        (
            momenti,
            new fvScalarMatrix
            (
                fvm::ddt(meanM)
              + meanMomentsAdvection_().divMoments()[momenti]
              - momentDiffusion(meanM)
            )
        );

        volUnivariateMoment& varM
        (
            meanMomentsVarianceQuadrature_.moments()[momenti]
        );

        meanMomentVarianceEqns.set
        (
            momenti,
            new fvScalarMatrix
            (
                fvm::ddt(varM)
              + meanMomentsVarianceAdvection_().divMoments()[momenti]
              - momentDiffusion(varM)
              ==
                envMixingModel_->K
                    (
                        meanM,
                        varM,
                        meanXi_
                    )
            )
        );
    }

    // Update moment environments to apply ODE solver
    calcEnvironmentMoments();

    // Solve source terms
    explicitMomentSource(envOneQuadrature_, 1);
    explicitMomentSource(envTwoQuadrature_, 2);

    // Update mixed moments
    calcMixedMoments();

    // Finish solving for moments
    forAll (meanMomentEqns, mEqni)
    {
        const volUnivariateMoment& meanM
        (
            meanMomentsQuadrature_.moments()[mEqni]
        );

        const volUnivariateMoment& varM
        (
            meanMomentsVarianceQuadrature_.moments()[mEqni]
        );

        meanMomentEqns[mEqni] -= fvc::ddt(meanM);
        meanMomentVarianceEqns[mEqni] -= fvc::ddt(varM);

        meanMomentVarianceEqns[mEqni].relax();
        meanMomentVarianceEqns[mEqni].solve();

        meanMomentEqns[mEqni].relax();
        meanMomentEqns[mEqni].solve();       
    }

    meanMomentsQuadrature_.updateQuadrature();
    meanMomentsVarianceQuadrature_.updateQuadrature();
}

// ************************************************************************* //
