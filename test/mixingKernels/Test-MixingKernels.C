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
    Test-MixingKernels

Description
    Checks the mixing kernels on the moments of a mixture fraction in a
    single cell of homogeneous turbulence, where their effect is known in
    closed form.

    The moments are those of a mixture fraction on [0, 1], so the zero-order
    moment is one and stays one. Both kernels conserve the mean, and both
    take the variance down exponentially at the rate the interaction by
    exchange with the mean gives it,

        d<xi'^2>/dt = -Cphi epsilon/k <xi'^2>,

    with Cphi = 2 the usual value. The Fokker-Planck kernel is built so that
    this holds whatever its coefficient Cmixing is, which only reshapes the
    distribution the variance belongs to: two of the cases run it with two
    values of Cmixing and have to decay at the same rate as the third, which
    runs the interaction by exchange with the mean.

    The population starts as almost all of it in the two pure streams, with
    one of its nodes at a mixture fraction of zero. That node used to be
    given no weight, as a size that has been driven to nothing is, which
    took the zero-order moment of the case from one to 0.2 in a step and
    left the kernels, written for a zero-order moment of one, destroying the
    mean.

    Run from the case directory.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "hexCellFvMesh.H"
#include "mixingModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    // The kernels take k and epsilon from the turbulence model, which is not
    // corrected here: the turbulence stays as the case gives it
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
        mixingProperties.subDict("turbulentMixingCoeffs")
       .subDict("mixingKernel");

    const scalar Cphi =
        kernelDict.getOrDefault<dimensionedScalar>
        (
            "Cphi",
            dimensionedScalar("Cphi", dimless, 2.0)
        ).value();

    const scalar epsilonByK =
        turbulence->epsilon()()[celli]/turbulence->k()()[celli];

    const volScalarField& m0 =
        mesh.lookupObject<volScalarField>("moment.0.mixing");

    const volScalarField& m1 =
        mesh.lookupObject<volScalarField>("moment.1.mixing");

    const volScalarField& m2 =
        mesh.lookupObject<volScalarField>("moment.2.mixing");

    const scalar mean0 = m1[celli]/m0[celli];
    const scalar variance0 = m2[celli]/m0[celli] - sqr(mean0);

    Info<< "\nMixing kernel " << kernelDict.get<word>("mixingKernel")
        << ", Cphi " << Cphi << ", epsilon/k " << epsilonByK << nl
        << "Mean " << mean0 << ", variance " << variance0 << nl << endl;


    // * * * * * * * * * * * * * * * * The steps * * * * * * * * * * * * * //

    const label nSteps = 20;

    // The zero-order moment and the mean only move by round-off, and the
    // variance is integrated by the solver to the tolerances of the case,
    // which take it to 1e-11 over the twenty steps. A kernel mixing at
    // twice the rate misses it by 0.9.
    const scalar conservationTolerance = 1.0e-10;
    const scalar decayTolerance = 1.0e-9;

    scalar worstM0 = 0;
    scalar worstMean = 0;
    scalar worstDecay = 0;

    for (label stepi = 0; stepi < nSteps; stepi++)
    {
        runTime++;
        mixing->solve();

        const scalar t = runTime.value();

        const scalar mean = m1[celli]/max(m0[celli], VSMALL);
        const scalar variance = m2[celli]/max(m0[celli], VSMALL) - sqr(mean);
        const scalar expected = variance0*Foam::exp(-Cphi*epsilonByK*t);

        const scalar errorM0 = mag(m0[celli] - 1.0);
        const scalar errorMean = mag(mean - mean0)/mean0;
        const scalar errorDecay = mag(variance - expected)/expected;

        Info<< "  t = " << t
            << ": m0 " << m0[celli]
            << ", mean " << mean
            << ", variance " << variance
            << " (expected " << expected
            << ", error " << errorDecay << ")" << endl;

        worstM0 = max(worstM0, errorM0);
        worstMean = max(worstMean, errorMean);
        worstDecay = max(worstDecay, errorDecay);
    }


    // * * * * * * * * * * * * * * * * The checks * * * * * * * * * * * * * //

    label nFailed = 0;

    Info<< nl;

    if (!(worstM0 <= conservationTolerance))
    {
        Info<< "FAILED: the zero-order moment moved from one by "
            << worstM0 << endl;
        nFailed++;
    }
    else
    {
        Info<< "OK: the zero-order moment stays one (" << worstM0 << ")"
            << endl;
    }

    if (!(worstMean <= conservationTolerance))
    {
        Info<< "FAILED: the mean moved by " << worstMean << endl;
        nFailed++;
    }
    else
    {
        Info<< "OK: the mean is conserved (" << worstMean << ")" << endl;
    }

    if (!(worstDecay <= decayTolerance))
    {
        Info<< "FAILED: the variance does not decay at Cphi epsilon/k, "
            << "largest relative error " << worstDecay << endl;
        nFailed++;
    }
    else
    {
        Info<< "OK: the variance decays at Cphi epsilon/k ("
            << worstDecay << ")" << endl;
    }

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
