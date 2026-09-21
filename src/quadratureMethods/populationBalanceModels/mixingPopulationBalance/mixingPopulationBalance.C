/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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
#include "supportType.H"

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
    odeType(phi.mesh(), dict),
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
        dict.lookupOrDefault<scalar>("minMixtureFractionVariance", 1.0e-4)
    ),
    minEnvironmentWeight_
    (
        dict.lookupOrDefault<scalar>("minEnvironmentWeight", 1.0e-6)
    ),
    p1_
    (
        mixingModel_().quadrature().nodes()[0].weight()
    ),
    xi1_
    (
        mixingModel_().quadrature().nodes()[0].abscissae()[0]
    ),
    p2_
    (
        mixingModel_().quadrature().nodes()[1].weight()
    ),
    xi2_
    (
        mixingModel_().quadrature().nodes()[1].abscissae()[0]
    ),
    meanXi_
    (
        mixingModel_().quadrature().moments()[1]
    ),
    meanMomentsQuadrature_
    (
        name + "MeanMoments",
        phi_.mesh(),
        List<supportType>(1, supportType::RPlus)
    ),
    meanMomentsVarianceQuadrature_
    (
        name + "MeanMomentsVariance",
        phi_.mesh(),
        List<supportType>(1, supportType::RPlus)
    ),
    meanMomentsAdvection_
    (
        univariateMomentAdvection::New
        (
            meanMomentsQuadrature_.subDict("momentAdvection"),
            meanMomentsQuadrature_,
            phi_,
            supportType::RPlus
        )
    ),
    meanMomentsVarianceAdvection_
    (
        univariateMomentAdvection::New
        (
            meanMomentsVarianceQuadrature_.subDict("momentAdvection"),
            meanMomentsVarianceQuadrature_,
            phi_,
            supportType::RPlus
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
    nucleation_(dict.lookup("nucleation")),
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
    growthModel_
    (
        Foam::populationBalanceSubModels::growthModel::New
        (
            dict.subDict("growthModel"),
            phi_.mesh()
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
::calcEnvironmentMoments()
{
    const volScalarMomentFieldSet& mXi
    (
        mixingModel_().quadrature().moments()
    );

    const volScalarMoment& xiMean_ = mXi[1];
    const volScalarMoment& xiMTwo_ = mXi[2];

    // Compute variance of the mixture fraction
    const volScalarField xiVariance(xiMTwo_ - sqr(xiMean_));

    // Difference between abscissae of the mixture fraction
    volScalarField xiDiff(xi1_ - xi2_);

    // The moments of the two environments are those that give back the
    // mean moments and the ones weighted by the mixture fraction,
    //
    //     p1 m1 + p2 m2 = <m>,  p1 xi1 m1 + p2 xi2 m2 = <xi m>,
    //
    // which takes the two environments to be apart: the solution divides
    // by p1, by p2 and by xi1 - xi2, with its sign. They are apart where
    // the variance of the mixture fraction, p1 p2 (xi1 - xi2)^2, is above
    // minMixtureFractionVariance, and each of them holds more than
    // minEnvironmentWeight of the fluid. Elsewhere the fluid of the cell is
    // mixed, to what those thresholds resolve, and both environments hold
    // the mean moments.
    //
    // The test used to be on xi1 - xi2 itself, which is negative wherever
    // the environments are apart, the nodes of the mixture fraction coming
    // in ascending order: no cell was ever split, and the model was one
    // of a single environment. An environment below minEnvironmentWeight
    // used to be given null moments where the other was split from it,
    // which would have taken what it held out of the mean.
    //
    // The solution is a population only where <xi m>/<m> lies between the
    // two mixture fractions, for every moment: outside, one environment
    // comes out negative, which transport, solving the two sets of moments
    // apart from the mixture fraction, can leave a cell with. Such a cell
    // is taken as mixed too: a negative environment used to be inverted to
    // nothing, and what it held, with its sign, went out of the mean. Which
    // cells hold two environments can be found from the fields of a time:
    // the two sets of moments and those of the mixture fraction say it.
    forAll(xiDiff, celli)
    {
        bool apart =
            xiVariance[celli] > minMixtureFractionVariance_
         && p1_[celli] > minEnvironmentWeight_
         && p2_[celli] > minEnvironmentWeight_;

        if (apart)
        {
            forAll(mEnvOne_, mi)
            {
                mEnvOne_[mi][celli] =
                    (
                        meanMomentsVariance_[mi][celli]
                      - meanMoments_[mi][celli]*xi2_[celli]
                    )/(p1_[celli]*xiDiff[celli]);

                mEnvTwo_[mi][celli] =
                    (
                        meanMoments_[mi][celli]*xi1_[celli]
                      - meanMomentsVariance_[mi][celli]
                    )/(p2_[celli]*xiDiff[celli]);

                // Round-off of a population that is all in one environment
                // is not a negative one
                const scalar small =
                    1e-12*mag(meanMoments_[mi][celli])
                   /min(p1_[celli], p2_[celli]);

                if
                (
                    mEnvOne_[mi][celli] < -small
                 || mEnvTwo_[mi][celli] < -small
                )
                {
                    apart = false;
                }
                else
                {
                    mEnvOne_[mi][celli] = max(mEnvOne_[mi][celli], scalar(0));
                    mEnvTwo_[mi][celli] = max(mEnvTwo_[mi][celli], scalar(0));
                }
            }
        }

        if (!apart)
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

        meanMomentsVariance_[mi] ==
            p1_*xi1_*mEnvOne_[mi] + p2_*xi2_*mEnvTwo_[mi];
    }
}


void
Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::updateCellMomentSource(const label)
{}


Foam::scalar
Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label environment
)
{
    scalar source = 0.0;
    if (aggregation_)
    {
        source +=
            aggregationKernel_->aggregationSource
            (
                momentOrder,
                celli,
                quadrature,
                environment
            );
    }
    if (breakup_)
    {
        source +=
            breakupKernel_->breakupSource
            (
                momentOrder,
                celli,
                quadrature
            );
    }
    if (growth_)
    {
        source +=
            growthModel_->phaseSpaceConvection
            (
                momentOrder,
                celli,
                quadrature
            );
    }
    if (nucleation_)
    {
        source += nucleationModel_->nucleationSource(momentOrder[0], celli);
    }

    return source;
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
    // The old-time moments exist before the sources overwrite the current
    // ones, so that the time derivative below is the one over the step
    forAll(meanMoments_, mi)
    {
        meanMoments_[mi].oldTime();
        meanMomentsVariance_[mi].oldTime();
    }

    // The population is shared out between the environments, changed by
    // its sources in each and put together again with the environments of
    // one time, the one the step starts from, which is also the time of the
    // moments shared out. The mixing model used to be solved first: the
    // moments weighted by the mixture fraction were then those of the old
    // environments and the mixture fractions those of the new, and a
    // population held by one stream alone, <xi m>/<m> being its mixture
    // fraction, fell outside the new ones by what the stream had mixed in
    // the step, which gave the other environment a negative population.
    Info << "Solving population balance" << endl;

    calcEnvironmentMoments();

    odeType::solve(envOneQuadrature_, 1);
    odeType::solve(envTwoQuadrature_, 2);

    calcMixedMoments();

    // Solve mixingModel
    Info << "Solving mixing model" << endl;

    mixingModel_().solve();

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
        volScalarMoment& meanM
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
              - diffusionModel_->momentDiff(meanM)
            )
        );

        volScalarMoment& varM
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
              - diffusionModel_->momentDiff(varM)
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

    // Finish solving for moments: what the sources changed enters as the
    // time derivative of the moments they overwrote
    forAll (meanMomentEqns, mEqni)
    {
        const volScalarMoment& meanM
        (
            meanMomentsQuadrature_.moments()[mEqni]
        );

        const volScalarMoment& varM
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


bool Foam::PDFTransportModels::populationBalanceModels::mixingPopulationBalance
::readIfModified()
{
    if (populationBalanceProperties_.modified())
    {
        odeType::read
        (
            populationBalanceProperties_.subDict(type() + "Coeffs")
        );
        return true;
    }

    return false;
}

// ************************************************************************* //
