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
    Copyright (C) 2019 Alberto Passalacqua
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
        mixingModel_().quadrature().nodes()[0].primaryWeight()
    ),
    xi1_
    (
        mixingModel_().quadrature().nodes()[0].primaryAbscissae()[0]
    ),
    p2_
    (
        mixingModel_().quadrature().nodes()[1].primaryWeight()
    ),
    xi2_
    (
        mixingModel_().quadrature().nodes()[1].primaryAbscissae()[0]
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

    // Compute moments in the environment
    forAll(xiDiff, celli)
    {
        // Null or very SMALL variance of the mixture fraction
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

    // Update moment environments to apply ODE solver
    calcEnvironmentMoments();

    // Solve source terms
    odeType::solve(envOneQuadrature_, 1);
    odeType::solve(envTwoQuadrature_, 2);

    // Update mixed moments
    calcMixedMoments();

    // Finish solving for moments
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
