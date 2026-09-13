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
    Test-RealizableOdeSolver

Description
    Checks that the realizable ODE solver integrates a source over every
    global time step, whatever the length of the steps and however small the
    moments are in their own units.

    A population of three sizes grows at a constant rate, which moves every
    size by Cg per unit time and leaves the moments known in closed form:
    with volume fractions w for weights and a length for abscissa, the
    number density of a size is n = w/x^3 and the moment of order k at time
    t is the sum of n (x + Cg t)^(k + 3) over the sizes.

    The case is driven through forty steps of four different lengths, so
    that the last substep of a step is as often as not a clipped one, and
    its weights are volume fractions of a dilute population, of order
    1e-10, so that what the solver takes as negligible has to be measured
    against the moments rather than against a fixed number. The solver used
    to declare a source null when it moved no moment by more than SMALL, an
    absolute number, and to carry a clipped substep of round-off from one
    step into the next as the step to start it with, on which every source
    was then null: a cell in that state integrated nothing from then on,
    and the moments here would stay where they started. It also began a
    global step with the local step it had ended the last one with, however
    long that was, so a step shorter than the last was integrated over the
    length of the last: the moments here would overshoot.

    Run from the testCase directory.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "populationBalanceModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // The global scope carries a second set of these, which makes a call
    // with plain scalars ambiguous
    using Foam::pow;

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // The population balance takes a flux, which nothing here transports
    // the moments with: the case lives in one cell and advects nothing
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("phi", dimVolume/dimTime, Zero)
    );

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
    const label nMoments = 6;

    // The population the case was built of, and the rate it grows at
    const scalarList w({2.0e-11, 5.0e-11, 3.0e-11});
    const scalarList x({1.0e-5, 2.0e-5, 4.0e-5});

    const scalar Cg =
        populationBalanceProperties.subDict("univariateCoeffs")
       .subDict("growthModel").get<scalar>("Cg");

    auto momentName = [](const label k)
    {
        return IOobject::groupName
        (
            IOobject::groupName("moment", Foam::name(k)), "populationBalance"
        );
    };

    auto expected = [&](const label k, const scalar t)
    {
        scalar m = 0;

        forAll(w, i)
        {
            m += w[i]/pow3(x[i])*pow(x[i] + Cg*t, k + 3);
        }

        return m;
    };


    // * * * * * * * * * * * * * * * * The steps * * * * * * * * * * * * * //

    // Four lengths of step, none a multiple of another, cycled through
    const scalarList deltaTs({1.0e-3, 1.7e-3, 0.6e-3, 1.3e-3});
    const label nSteps = 40;

    Info<< "\nGrowing the population over " << nSteps
        << " steps of four lengths" << endl;

    for (label stepi = 0; stepi < nSteps; stepi++)
    {
        runTime.setDeltaT(deltaTs[stepi % deltaTs.size()]);
        runTime++;
        populationBalance->solve();
    }

    const scalar t = runTime.value();


    // * * * * * * * * * * * * * * * * The checks * * * * * * * * * * * * * //

    Info<< "\nThe moments after " << t << " s" << endl;

    // The moments are polynomials in time of degree up to eight, which the
    // third-order scheme integrates to a truncation error that a step of
    // a thousandth of the size of a particle leaves far below this
    const scalar tolerance = 1.0e-7;

    label nFailed = 0;

    for (label k = 0; k < nMoments; k++)
    {
        const scalar computed =
            mesh.lookupObject<volScalarField>(momentName(k))[celli];

        const scalar reference = expected(k, t);
        const scalar start = expected(k, 0);

        // The growth the population has undergone, as the fraction of the
        // moment it changed by: what a frozen cell would miss entirely
        const scalar growth = mag(reference - start)/mag(reference);
        const scalar error = mag(computed - reference)/mag(reference);

        Info<< "  moment " << k << ": " << computed
            << ", expected " << reference
            << " (grew by " << growth << ", error " << error << ")" << endl;

        if (!(error <= tolerance))
        {
            nFailed++;
        }
    }

    if (nFailed > 0)
    {
        FatalErrorInFunction
            << nFailed << " of " << nMoments << " moments were not "
            << "integrated to the growth the population underwent." << nl
            << exit(FatalError);
    }

    Info<< "\nEvery moment was integrated over every step." << nl << endl;

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
