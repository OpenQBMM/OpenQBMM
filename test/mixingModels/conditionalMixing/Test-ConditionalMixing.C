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
    Test-ConditionalMixing

Description
    Test the mixing of a mixture fraction and of a composition conditioned
    on it, in one cell of homogeneous turbulence, where what mixing does is
    known.

    The distribution is three nodes of the mixture fraction, each carrying
    one composition. IEM relaxes every value to its mean at Cphi/2
    epsilon/k, so over the steps:

    - the zero-order moment stays one and the means of the mixture fraction
      and of the composition are conserved, since mixing moves neither;
    - the variance of the mixture fraction decays at Cphi epsilon/k;
    - the covariance of the mixture fraction and the composition decays at
      the same rate, which is what the source of the moment of orders
      (1, 1) comes to;
    - every conditional mean stays within [0, 1].

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "hexCellFvMesh.H"
#include "mixingModel.H"
#include "quadratureApproximations.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nChecks = 0;
label nFailed = 0;

void check(const word& name, const scalar error, const scalar tolerance)
{
    nChecks++;

    if (error <= tolerance)
    {
        Info<< "    passed: " << name << " (" << error << ")" << nl;
    }
    else
    {
        nFailed++;
        Info<< "    FAILED: " << name << " (" << error << ", tolerance "
            << tolerance << ")" << nl;
    }
}


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

    // The kernel takes k and epsilon from the turbulence model, which is
    // not corrected here: the turbulence stays as the case gives it
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New(rho, U, phi, thermo)
    );

    turbulence->validate();

    IOdictionary mixingProperties
    (
        IOobject
        (
            "mixingProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<mixingModel> mixing
    (
        mixingModel::New("mixing", mixingProperties, phi)
    );

    const label celli = 0;

    const dictionary& kernelDict =
        mixingProperties.subDict("conditionalMixingCoeffs")
       .subDict("mixingKernel");

    const scalar Cphi =
        kernelDict.getOrDefault<dimensionedScalar>
        (
            "Cphi",
            dimensionedScalar("Cphi", dimless, 2.0)
        ).value();

    const scalar epsilonByK =
        turbulence->epsilon()()[celli]/turbulence->k()()[celli];

    auto moment = [&mesh, celli](const word& orders)
    {
        return
            mesh.lookupObject<volScalarField>("moment." + orders + ".mixing")
            [celli];
    };

    const scalar m0 = moment("0_0");
    const scalar mean0 = moment("1_0")/m0;
    const scalar meanY0 = moment("0_1")/m0;
    const scalar variance0 = moment("2_0")/m0 - sqr(mean0);
    const scalar covariance0 = moment("1_1")/m0 - mean0*meanY0;

    Info<< "\nMixing " << kernelDict.get<word>("mixingKernel")
        << ", Cphi " << Cphi << ", epsilon/k " << epsilonByK << nl
        << "Mixture fraction " << mean0 << " of variance " << variance0
        << ", composition " << meanY0 << " of covariance " << covariance0
        << nl << endl;

    const velocityQuadratureApproximation& quadrature =
        mesh.lookupObject<velocityQuadratureApproximation>
        (
            "quadratureProperties.mixing"
        );


    // * * * * * * * * * * * * * * * * The steps * * * * * * * * * * * * * //

    const label nSteps = 20;

    scalar worstM0 = 0;
    scalar worstMean = 0;
    scalar worstMeanY = 0;
    scalar worstVariance = 0;
    scalar worstCovariance = 0;
    scalar worstBound = 0;

    for (label stepi = 0; stepi < nSteps; stepi++)
    {
        runTime++;
        mixing->solve();

        const scalar t = runTime.value();
        const scalar decay = Foam::exp(-Cphi*epsilonByK*t);

        const scalar m = moment("0_0");
        const scalar mean = moment("1_0")/m;
        const scalar meanY = moment("0_1")/m;
        const scalar variance = moment("2_0")/m - sqr(mean);
        const scalar covariance = moment("1_1")/m - mean*meanY;

        worstM0 = max(worstM0, mag(m - 1));
        worstMean = max(worstMean, mag(mean - mean0)/mean0);
        worstMeanY = max(worstMeanY, mag(meanY - meanY0)/meanY0);

        worstVariance =
            max(worstVariance, mag(variance - variance0*decay)/(variance0*decay));

        worstCovariance =
            max
            (
                worstCovariance,
                mag(covariance - covariance0*decay)/mag(covariance0*decay)
            );

        forAll(quadrature.nodes(), nodei)
        {
            const scalar y = quadrature.nodes()[nodei].abscissae()[1][celli];

            worstBound = max(worstBound, max(-y, y - 1));
        }

        Info<< "  t = " << t << ": variance " << variance
            << " (expected " << variance0*decay << "), covariance "
            << covariance << " (expected " << covariance0*decay << ")"
            << endl;
    }


    // * * * * * * * * * * * * * * * * The checks * * * * * * * * * * * * * //

    Info<< nl;

    // The moments that do not move are conserved to round-off; the two
    // that decay are integrated by the ODE solver to the tolerances of the
    // case, which take them to 1e-9 over the twenty steps
    check("the zero-order moment stays one", worstM0, 1e-10);
    check("the mean of the mixture fraction is conserved", worstMean, 1e-10);
    check("the mean of the composition is conserved", worstMeanY, 1e-10);

    check
    (
        "the variance decays at Cphi epsilon/k",
        worstVariance,
        1e-8
    );

    check
    (
        "the covariance decays at Cphi epsilon/k",
        worstCovariance,
        1e-8
    );

    check("every conditional mean stays within [0, 1]", worstBound, 0);

    if (nFailed > 0)
    {
        FatalErrorInFunction
            << nFailed << " of " << nChecks << " checks failed." << nl
            << exit(FatalError);
    }

    Info<< nl << "All " << nChecks << " checks passed" << nl << endl;

    return 0;
}


// ************************************************************************* //
