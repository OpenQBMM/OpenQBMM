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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-SizeVelocityNucleation

Description
    Checks the source the nucleation of particles adds to the moments of a
    size-velocity population balance.

    The nuclei form out of the continuous phase and, at the size they
    appear with, follow it without slip, so they are born with its
    velocity: in velocity space the source is a delta at Uc. The moment of
    order (k, i, j) therefore gains

        J*x^k*Uc_x^i*Uc_y^j

    per unit time, where J is the rate of the nucleation model and x the
    size of a nucleus.

    The weights of this case are volume fractions and its size coordinate
    is a diameter, so its moments carry the volume of the particles rather
    than their number, and the order of the size is three above the one the
    moment is named by: the nuclei are added in volume as well, as
    aggregation, breakup and growth already add theirs.

    The case runs one step of the population balance with nucleation as its
    only source and the explicit Euler of the realizable ODE solver, so
    that what a moment gains over the step is its source times the length
    of the step. The rate and the size of the nuclei are read here from the
    dictionary of the model, so that the closure is checked against them
    rather than against itself.

    Run from the testCase directory.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "populationBalanceModel.H"
#include "fundamentalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nChecks = 0;
label nFailed = 0;


//- Compare a computed value with the one expected, relative to a scale
void compare
(
    const scalar computed,
    const scalar expected,
    const scalar scale,
    const string& what,
    const scalar tolerance = 1.0e-8
)
{
    nChecks++;

    const scalar error = mag(computed - expected)/max(mag(scale), VSMALL);

    // A value that is not a number fails the first comparison and an
    // infinite one the second
    if (computed == computed && mag(computed) <= VGREAT && error <= tolerance)
    {
        Info<< "  OK: \"" << what.c_str() << "\" (relative error "
            << error << ")" << nl;

        return;
    }

    nFailed++;

    Info<< "  FAILED: \"" << what.c_str() << "\"" << nl
        << "      computed " << computed << ", expected " << expected
        << ", relative error " << error << nl;
}


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

    // The velocity of the gas, which the nuclei are born with. The
    // population balance looks it up rather than being given it.
    volVectorField Uc
    (
        IOobject
        (
            "U.gas",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // The population balance takes a flux, which nothing here transports
    // the moments with: the case advects nothing and lives in one cell
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

    // The orders of the moments the case carries, as its
    // quadratureProperties lists them
    const labelListList momentOrders
    (
        IOdictionary
        (
            IOobject
            (
                "quadratureProperties.populationBalance",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).lookup("moments")
    );

    auto momentName = [](const labelList& order)
    {
        word w;

        forAll(order, dimi)
        {
            w += Foam::name(order[dimi]);
        }

        return IOobject::groupName
        (
            IOobject::groupName("moment", w), "populationBalance"
        );
    };

    scalarList before(momentOrders.size(), Zero);

    forAll(momentOrders, mi)
    {
        before[mi] =
            mesh.lookupObject<volScalarField>
            (
                momentName(momentOrders[mi])
            )[celli];
    }


    // * * * * * * * * * * * * * * * * The step * * * * * * * * * * * * * * //

    runTime++;

    populationBalance->solve();

    const scalar dt = runTime.deltaTValue();

    scalarList after(momentOrders.size(), Zero);

    forAll(momentOrders, mi)
    {
        after[mi] =
            mesh.lookupObject<volScalarField>
            (
                momentName(momentOrders[mi])
            )[celli];
    }


    // * * * * * * * * * * * * * * What is expected * * * * * * * * * * * * //

    // The rate of the nucleation model and the size of a nucleus, read
    // here from the dictionary of the case rather than from the model
    const dictionary& nucleationDict =
        populationBalanceProperties.subDict("sizeVelocityCoeffs")
       .subDict("nucleationModel");

    const scalar J =
        dimensionedScalar("nucleationRate", nucleationDict).value();

    const scalar xNucleation =
        dimensionedScalar("nucleationSize", nucleationDict).value();

    Info<< "\nNuclei of " << xNucleation << " m, "
        << J << " of them per cubic metre and second" << nl;


    // * * * * * * * * * * * * * * * * The checks * * * * * * * * * * * * * //

    Info<< "\nThe source of nucleation" << endl;

    forAll(momentOrders, mi)
    {
        const labelList& order = momentOrders[mi];

        // The weights are volume fractions and the size is a diameter,
        // so the order of the size is three above the one of the moment
        scalar expected = dt*J*pow(xNucleation, order[0] + 3);

        for (label cmpt = 0; cmpt < 2; cmpt++)
        {
            expected *= pow(Uc[celli][cmpt], order[cmpt + 1]);
        }

        compare
        (
            after[mi] - before[mi],
            expected,
            expected,
            "the moment of order " + Foam::name(order[0])
          + " " + Foam::name(order[1]) + " " + Foam::name(order[2])
          + " gains J*x^k*Uc^(i,j)"
        );
    }

    // The two properties the closure rests on, stated on their own so that
    // a failure says which of them broke
    compare
    (
        after[0] - before[0],
        dt*J*pow(xNucleation, 3),
        dt*J*pow(xNucleation, 3),
        "nucleation adds volume where the weights are volume fractions"
    );

    {
        // The velocity the nuclei are born with, recovered from the two
        // moments of order one in velocity
        const label iUx = momentOrders.find(labelList({0, 1, 0}));
        const label iUy = momentOrders.find(labelList({0, 0, 1}));

        compare
        (
            (after[iUx] - before[iUx])/(after[0] - before[0]),
            Uc[celli].x(),
            Uc[celli].x(),
            "the nuclei are born with the velocity of the continuous phase, "
            "first component"
        );

        compare
        (
            (after[iUy] - before[iUy])/(after[0] - before[0]),
            Uc[celli].y(),
            Uc[celli].y(),
            "the nuclei are born with the velocity of the continuous phase, "
            "second component"
        );
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl;

    if (nFailed > 0)
    {
        FatalErrorInFunction
            << nFailed << " of " << nChecks << " checks failed."
            << exit(FatalError);
    }

    Info<< nChecks << " checks passed." << nl << endl;

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
