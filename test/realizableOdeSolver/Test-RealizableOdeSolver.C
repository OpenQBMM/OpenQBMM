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
    moments are in their own units, and that it subdivides a step its
    tolerances ask it to.

    A population of three sizes moves at the rate the case selects, which
    leaves the moments known in closed form: with volume fractions w for
    weights and a length for abscissa, the number density of a size is
    n = w/x^3 and the moment of order k at time t is the sum of
    n x(t)^(k + 3) over the sizes, with x(t) the size the model gives.

    The cases are run from their own directories:

    - testCase, a constant growth rate, which moves every size by Cg per
      unit time. The moments are polynomials in time, integrated exactly by
      the third-order scheme, so the answer is right to round-off unless
      the solver drops part of a step.

    - nucleationCase, particles of one size formed at a constant rate, and
      no growth. The source of every moment is the same at every stage of a
      step, so the error the solver estimates is exactly zero, on which the
      factor it grows the next step by used to be a division by zero.

    - stiffCase, the non-linear evaporation of the same population over
      steps ten times longer, which takes the smallest size down to a
      quarter of what it started as. The rate diverges as a droplet
      vanishes, so a step that is a fair one at the start is far too long
      by the end, and the solver has to see that and subdivide. Its
      tolerances are far above every moment of the case as absolute
      numbers: measured against the scale of each moment they are met at
      4e-8, and taken as plain numbers, as they were, the controller sees
      no error at all, integrates every step whole and ends 1e-4 away.

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
    using Foam::sqrt;

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

    // The population the case was built of
    const scalarList w({2.0e-11, 5.0e-11, 3.0e-11});
    const scalarList x({1.0e-5, 2.0e-5, 4.0e-5});

    const dictionary& coeffsDict =
        populationBalanceProperties.subDict("univariateCoeffs");

    // The rate the sizes grow at, where the case grows them
    const bool growth = coeffsDict.get<Switch>("growth");

    word growthModel;
    scalar Cg = 0;

    if (growth)
    {
        const dictionary& growthDict = coeffsDict.subDict("growthModel");

        growthModel = growthDict.get<word>("growthModel");
        Cg = growthDict.get<scalar>("Cg");
    }

    // The rate particles form at, and the size they form with, where the
    // case forms them
    const bool nucleation = coeffsDict.get<Switch>("nucleation");

    scalar J = 0;
    scalar xNuclei = 0;

    if (nucleation)
    {
        const dictionary& nucleationDict =
            coeffsDict.subDict("nucleationModel");

        J = dimensionedScalar("nucleationRate", nucleationDict).value();
        xNuclei = dimensionedScalar("nucleationSize", nucleationDict).value();
    }

    // The size of a droplet at a time. A constant rate moves it by Cg t.
    // The non-linear evaporation of a length coordinate is the d-square
    // law, dx/dt = -Cg/(2 pi x), which takes x^2 down by Cg t / pi and
    // whose rate diverges as a droplet vanishes: a step that is a fair
    // one at the size a droplet starts at is far too long by the time it
    // has evaporated, which is what the controller is there to see.
    auto size = [&](const scalar x0, const scalar t)
    {
        if (!growth)
        {
            return x0;
        }

        return
            growthModel == "constant"
          ? x0 + Cg*t
          : sqrt
            (
                max(sqr(x0) - Cg*t/constant::mathematical::pi, scalar(0))
            );
    };

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
            m += w[i]/pow3(x[i])*pow(size(x[i], t), k + 3);
        }

        // The weights are volume fractions, so the particles that formed
        // add their number times their size to the power k + 3
        m += J*t*pow(xNuclei, k + 3);

        return m;
    };


    // * * * * * * * * * * * * * * * * The steps * * * * * * * * * * * * * //

    // Four lengths of step, none a multiple of another, cycled through, in
    // units of the step the case gives: a case whose source is stiff asks
    // for a longer one, which the solver is then to subdivide itself
    const scalar deltaT0 = runTime.deltaTValue();

    const scalarList deltaTs
    ({
        1.0*deltaT0, 1.7*deltaT0, 0.6*deltaT0, 1.3*deltaT0
    });

    const label nSteps = 40;

    Info<< "\nMoving the population over " << nSteps
        << " steps of four lengths, the longest " << deltaTs[1] << " s"
        << endl;

    for (label stepi = 0; stepi < nSteps; stepi++)
    {
        runTime.setDeltaT(deltaTs[stepi % deltaTs.size()]);
        runTime++;
        populationBalance->solve();
    }

    const scalar t = runTime.value();


    // * * * * * * * * * * * * * * * * The checks * * * * * * * * * * * * * //

    Info<< "\nThe moments after " << t << " s" << endl;

    // Constant growth leaves the moments polynomials in time of degree up
    // to eight, which the third-order scheme integrates to a truncation
    // error far below this; evaporation is met at this accuracy only by a
    // solver that subdivides the step it is given
    const scalar tolerance = 1.0e-7;

    label nFailed = 0;

    for (label k = 0; k < nMoments; k++)
    {
        const scalar computed =
            mesh.lookupObject<volScalarField>(momentName(k))[celli];

        const scalar reference = expected(k, t);
        const scalar start = expected(k, 0);

        // The change the population has undergone, as the fraction of the
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
