/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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
    along with OpenQBMM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-MixingPopulationBalance

Description
    Checks the population balance of two environments of a mixture fraction,
    mixingPbe, in a single cell of homogeneous turbulence, on a population
    whose evolution is known in closed form.

    The cell holds two streams, 0.7 of its fluid at a mixture fraction of
    zero and 0.3 at one, mixed by the interaction by exchange with the mean.
    The stream at one carries a population, of moments mu_k, and the other
    none, and nothing acts on the population but mixing. The population of
    an environment is then mu_k times its mixture fraction at every time,
    since both are mixed in the same way, which gives the two moment sets
    the model transports:

        <m_k> = mu_k <xi>,  <xi m_k> = mu_k <xi^2>,

    the first constant and the second following the second moment of the
    mixture fraction, that is its variance, which decays at Cphi epsilon/k.

    Two things of the model are held by it. The environments have to be
    told apart: where they are not, both take the mean moments and
    <xi m_k> is set to <xi><m_k>, that is mu_k <xi>^2, at the first step.
    They never were, the test that tells them apart having been on the
    difference of the two mixture fractions with its sign, which is
    negative. And the source that mixes <xi m_k> has to work at the rate of
    the mixing of the mixture fraction itself, which it did at twice that.

    <xi m_k> is solved implicitly with the Euler scheme, which leaves an
    error of first order in the time step on it: a part in a thousand here.

    Run from the case directory.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "hexCellFvMesh.H"
#include "populationBalanceModel.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"

    Foam::simplifiedMeshes::hexCellFvMesh mesh(runTime);

    autoPtr<psiThermo> pThermo(psiThermo::New(mesh));
    psiThermo& thermo = pThermo();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.rho()
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, Zero),
        thermo.p().boundaryField().types()
    );

    #include "createPhi.H"

    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New(rho, U, phi, thermo)
    );

    turbulence->validate();

    IOdictionary populationBalanceProperties
    (
        IOobject
        (
            "populationBalanceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<populationBalanceModel> populationBalance
    (
        populationBalanceModel::New
        (
            "populationBalance", populationBalanceProperties, phi
        )
    );

    const label celli = 0;
    const label nMoments = 4;

    const scalar Cphi =
        populationBalanceProperties.subDict("mixingPbeCoeffs")
       .subDict("environmentMixingModel")
       .get<dimensionedScalar>("Cphi").value();

    const scalar epsilonByK =
        turbulence->epsilon()()[celli]/turbulence->k()()[celli];

    UPtrList<const volScalarField> xi(3);
    UPtrList<const volScalarField> mean(nMoments);
    UPtrList<const volScalarField> weighted(nMoments);

    forAll(xi, k)
    {
        xi.set
        (
            k,
            &mesh.lookupObject<volScalarField>
            (
                "moment." + Foam::name(k) + ".mixing"
            )
        );
    }

    forAll(mean, k)
    {
        mean.set
        (
            k,
            &mesh.lookupObject<volScalarField>
            (
                "moment." + Foam::name(k) + ".populationBalanceMeanMoments"
            )
        );

        weighted.set
        (
            k,
            &mesh.lookupObject<volScalarField>
            (
                "moment." + Foam::name(k)
              + ".populationBalanceMeanMomentsVariance"
            )
        );
    }

    // The population of the stream at one, and the state of the mixture
    // fraction the run starts from
    const scalar xiMean0 = xi[1][celli];
    const scalar xiVariance0 = xi[2][celli] - sqr(xiMean0);

    scalarList mu(nMoments);

    forAll(mu, k)
    {
        mu[k] = mean[k][celli]/xiMean0;
    }

    Info<< "\nTwo streams of mean mixture fraction " << xiMean0
        << " and variance " << xiVariance0 << ", Cphi epsilon/k "
        << Cphi*epsilonByK << nl
        << "Population of the stream at one: " << mu << nl << endl;

    scalar worstMean = 0;
    scalar worstWeighted = 0;
    scalar worstXi = 0;

    while (runTime.run())
    {
        runTime++;
        populationBalance->solve();

        const scalar t = runTime.value();

        // The second moment of the mixture fraction, exact and computed
        const scalar xiTwoExact =
            sqr(xiMean0) + xiVariance0*Foam::exp(-Cphi*epsilonByK*t);

        worstXi = max(worstXi, mag(xi[2][celli] - xiTwoExact)/xiTwoExact);

        forAll(mu, k)
        {
            worstMean =
                max
                (
                    worstMean,
                    mag(mean[k][celli] - mu[k]*xiMean0)/(mu[k]*xiMean0)
                );

            worstWeighted =
                max
                (
                    worstWeighted,
                    mag(weighted[k][celli] - mu[k]*xiTwoExact)
                   /(mu[k]*xiTwoExact)
                );
        }

        if (runTime.timeIndex() % 100 == 0)
        {
            Info<< "  t = " << t << ": <xi^2> " << xi[2][celli]
                << " (exact " << xiTwoExact << "), <xi m_0>/mu_0 "
                << weighted[0][celli]/mu[0] << endl;
        }
    }

    label nFailed = 0;

    Info<< nl;

    auto check = [&nFailed]
    (
        const string& what,
        const scalar error,
        const scalar tolerance
    )
    {
        if (!(error <= tolerance))
        {
            Info<< "FAILED: " << what.c_str() << ", largest relative error "
                << error << " (tolerance " << tolerance << ")" << endl;

            nFailed++;
        }
        else
        {
            Info<< "OK: " << what.c_str() << " (" << error << ")" << endl;
        }
    };

    check
    (
        "the second moment of the mixture fraction decays at Cphi epsilon/k",
        worstXi,
        1e-8
    );

    check
    (
        "the mean moments of a population that only mixes are constant",
        worstMean,
        1e-10
    );

    check
    (
        "the moments weighted by the mixture fraction follow its second "
        "moment",
        worstWeighted,
        2e-3
    );

    if (nFailed > 0)
    {
        FatalErrorInFunction
            << nFailed << " of 3 checks failed." << nl
            << exit(FatalError);
    }

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
