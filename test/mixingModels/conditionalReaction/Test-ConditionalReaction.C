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
    Test-ConditionalReaction

Description
    Test the chemistry of OpenFOAM solved in every environment of the
    mixture fraction, in one cell, where what the reaction does is known.

    The distribution is three nodes of the mixture fraction, each carrying
    the mass fractions of two species A and B, with Y_A = xi and Y_B = 1 -
    xi at every node, and the reaction A = B proceeds at a constant rate k
    whatever the temperature. Since a reaction moves neither the weights
    nor the mixture fraction, every moment <xi^k Y_A> decays as exp(-k t)
    and <xi^k Y_B> takes what it loses. The case says what else is on:
    - with no mixing, the moments of the mixture fraction do not change
      and every <xi^k Y_A>, up to the second order, decays exactly;
    - with IEM, the mean of Y_A still decays exactly, since mixing
      conserves it, and the zeroth moment and the mean of the mixture
      fraction are conserved;
    - with the enthalpy linear in the mixture fraction, the temperature of
      the cell is the mean of the temperatures of the environments, each
      the one at which its composition has the enthalpy of its mixture
      fraction, and it rises, since A carries a heat of formation;
    - with isothermal environments, the temperature of the cell does not
      change;
    - with the enthalpy conditioned on the mixture fraction, a coordinate
      of the distribution that starts linear in it between the streams,
      the enthalpy stays on that line while the environments mix, since
      IEM keeps a linear relation linear, and the temperatures are those
      the linear closure gives. The moments of the enthalpy stay on it to
      the tolerance of the solver of the moment equations. The enthalpy of
      the nodes carries that error amplified by the solve for the
      conditional means, more and more as mixing draws the nodes together:
      a hundred times and more by the end of the test, which the tolerance
      on the nodes allows for.
    In every case the mass fractions of every node stay within [0, 1] and
    sum to one, the species fields of the thermo hold the mean composition,
    and the temperature field of every environment holds the temperature
    the closure gives it.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "psiReactionThermo.H"
#include "basicSpecieMixture.H"
#include "turbulentFluidThermoModel.H"
#include "hexCellFvMesh.H"
#include "mixingModel.H"
#include "quadratureApproximations.H"
#include "IFstream.H"

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

    const basicSpecieMixture& composition =
        refCast<psiReactionThermo>(thermo).composition();

    const label A = composition.species().find("A");
    const label B = composition.species().find("B");

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

    const dictionary& coeffs =
        mixingProperties.subDict("conditionalMixingCoeffs");

    const bool mixingOn =
        coeffs.subDict("mixingKernel").get<word>("mixingKernel") != "none";

    const dictionary& reactionDict = coeffs.subDict("reaction");
    const word temperature(reactionDict.get<word>("temperature"));
    const bool conditionalEnthalpy = temperature == "conditionalEnthalpy";
    const bool linearEnthalpy =
        temperature == "linearEnthalpy" || conditionalEnthalpy;

    // The moments have one index more where the enthalpy is carried
    const word tail(conditionalEnthalpy ? "_0" : "");

    // The rate of A = B, from the reactions of the case
    IFstream reactionsFile(runTime.constant()/"reactions");
    const dictionary reactions(reactionsFile);
    const scalar k =
        reactions.subDict("reactions").subDict("AtoB").get<scalar>("A");

    // The streams the enthalpy is linear between
    scalar T0 = 0;
    scalar T1 = 0;
    scalarField Y0(2, Zero);
    scalarField Y1(2, Zero);

    if (linearEnthalpy)
    {
        const dictionary& streams = reactionDict.subDict("streams");
        T0 = streams.subDict("xi0").get<scalar>("T");
        T1 = streams.subDict("xi1").get<scalar>("T");
        const dictionary& YDict0 = streams.subDict("xi0").subDict("Y");
        const dictionary& YDict1 = streams.subDict("xi1").subDict("Y");
        Y0[A] = YDict0.getOrDefault<scalar>("A", 0);
        Y0[B] = YDict0.getOrDefault<scalar>("B", 0);
        Y1[A] = YDict1.getOrDefault<scalar>("A", 0);
        Y1[B] = YDict1.getOrDefault<scalar>("B", 0);
    }

    auto moment = [&mesh, &tail, celli](const word& orders)
    {
        return
            mesh.lookupObject<volScalarField>
            (
                "moment." + orders + tail + ".mixing"
            )[celli];
    };

    const velocityQuadratureApproximation& quadrature =
        mesh.lookupObject<velocityQuadratureApproximation>
        (
            "quadratureProperties.mixing"
        );

    // Absolute enthalpy of a composition at a temperature
    auto enthalpy = [&composition, &thermo, celli]
    (
        const scalarField& Y,
        const scalar T
    )
    {
        const scalar p = thermo.p()[celli];
        scalar h = 0;

        forAll(Y, j)
        {
            h += Y[j]*composition.Ha(j, p, T);
        }

        return h;
    };

    // The temperature at which a composition has an enthalpy
    auto temperatureOf = [&composition, &thermo, &enthalpy, celli]
    (
        const scalarField& Y,
        const scalar h,
        scalar T
    )
    {
        const scalar p = thermo.p()[celli];

        for (label iter = 0; iter < 100; iter++)
        {
            scalar Cp = 0;

            forAll(Y, j)
            {
                Cp += Y[j]*composition.Cp(j, p, T);
            }

            const scalar dT = (enthalpy(Y, T) - h)/Cp;
            T -= dT;

            if (mag(dT) < 1e-12*T)
            {
                break;
            }
        }

        return T;
    };

    const wordList YAOrders({"0_1_0", "1_1_0", "2_1_0"});
    const wordList YBOrders({"0_0_1", "1_0_1", "2_0_1"});
    const wordList xiOrders
    ({
        "0_0_0", "1_0_0", "2_0_0", "3_0_0", "4_0_0", "5_0_0"
    });

    scalarList YA0(3);
    scalarList xi0(6);

    forAll(YAOrders, i)
    {
        YA0[i] = moment(YAOrders[i]);
    }

    forAll(xiOrders, i)
    {
        xi0[i] = moment(xiOrders[i]);
    }

    const scalar Tinitial = thermo.T()[celli];

    Info<< "\nA = B at " << k << " 1/s, mixing "
        << (mixingOn ? "on" : "off") << ", temperature " << temperature
        << nl << "<Y_A> " << YA0[0] << ", <xi Y_A> " << YA0[1]
        << ", <xi^2 Y_A> " << YA0[2] << ", T " << Tinitial << " K"
        << nl << endl;


    // * * * * * * * * * * * * * * * * The steps * * * * * * * * * * * * * //

    const label nSteps = 20;

    scalar worstDecay = 0;
    scalar worstMeanDecay = 0;
    scalar worstSum = 0;
    scalar worstXi = 0;
    scalar worstM0 = 0;
    scalar worstMean = 0;
    scalar worstBound = 0;
    scalar worstMeanY = 0;
    scalar worstT = 0;
    scalar worstNodeT = 0;
    scalar worstLine = 0;
    scalar worstMomentLine = 0;
    scalar Tfinal = Tinitial;

    for (label stepi = 0; stepi < nSteps; stepi++)
    {
        runTime++;
        // The model gives the thermo the mean state and corrects it
        mixing->solve();

        const scalar t = runTime.value();
        const scalar decay = Foam::exp(-k*t);

        // The moments of A decay, and B takes what A loses
        forAll(YAOrders, i)
        {
            const scalar YA = moment(YAOrders[i]);
            const scalar YB = moment(YBOrders[i]);
            const scalar xi = moment(xiOrders[i]);

            worstSum = max(worstSum, mag(YA + YB - xi));

            const scalar error = mag(YA - decay*YA0[i])/YA0[i];

            if (i == 0)
            {
                worstMeanDecay = max(worstMeanDecay, error);
            }
            else if (!mixingOn)
            {
                worstDecay = max(worstDecay, error);
            }
        }

        // The mixture fraction
        forAll(xiOrders, i)
        {
            const scalar error = mag(moment(xiOrders[i]) - xi0[i]);

            if (i == 0)
            {
                worstM0 = max(worstM0, error);
            }
            else if (i == 1)
            {
                worstMean = max(worstMean, error);
            }
            else if (!mixingOn)
            {
                worstXi = max(worstXi, error);
            }
        }

        // The moments of the enthalpy, linear in the mixture fraction
        if (conditionalEnthalpy)
        {
            const scalar h0 = enthalpy(Y0, T0);
            const scalar h1 = enthalpy(Y1, T1);
            const wordList hOrders({"0_0_0_1", "1_0_0_1", "2_0_0_1"});

            forAll(hOrders, i)
            {
                const scalar mh =
                    mesh.lookupObject<volScalarField>
                    (
                        "moment." + hOrders[i] + ".mixing"
                    )[celli];

                const scalar exact =
                    h0*moment(xiOrders[i])
                  + (h1 - h0)*moment(xiOrders[i + 1]);

                worstMomentLine =
                    max(worstMomentLine, mag(mh - exact)/mag(h1 - h0));
            }
        }

        // The nodes: bounds, and the mean state the thermo received
        const auto& nodes = quadrature.nodes();

        scalar wSum = 0;
        scalar wT = 0;
        scalarField meanY(2, Zero);

        forAll(nodes, nodei)
        {
            const scalar w = nodes[nodei].weight()[celli];
            const scalar xi = nodes[nodei].abscissae()[0][celli];

            scalarField Y(2);
            Y[A] = nodes[nodei].abscissae()[1 + A][celli];
            Y[B] = nodes[nodei].abscissae()[1 + B][celli];

            forAll(Y, j)
            {
                worstBound = max(worstBound, max(-Y[j], Y[j] - 1));
            }

            worstSum = max(worstSum, mag(sum(Y) - 1));

            scalar T = Tinitial;

            if (linearEnthalpy)
            {
                const scalar h0 = enthalpy(Y0, T0);
                const scalar h1 = enthalpy(Y1, T1);
                scalar h = h0 + xi*(h1 - h0);

                if (conditionalEnthalpy)
                {
                    // The enthalpy the node carries, which has to be the
                    // one of the line
                    const scalar hNode = nodes[nodei].abscissae()[3][celli];

                    worstLine =
                        max(worstLine, mag(hNode - h)/mag(h1 - h0));

                    h = hNode;
                }

                T = temperatureOf(Y, h, Tinitial);
            }

            // The temperature field of the environment
            const volScalarField& Tnode =
                mesh.lookupObject<volScalarField>
                (
                    IOobject::groupName("T", nodes[nodei].name())
                );

            worstNodeT = max(worstNodeT, mag(Tnode[celli] - T)/T);

            wSum += w;
            wT += w*T;
            meanY += w*Y;
        }

        meanY /= wSum;
        const scalar Tmean = wT/wSum;

        worstMeanY =
            max
            (
                worstMeanY,
                max
                (
                    mag(composition.Y(A)[celli] - meanY[A]),
                    mag(composition.Y(B)[celli] - meanY[B])
                )
            );

        worstT = max(worstT, mag(thermo.T()[celli] - Tmean)/Tmean);
        Tfinal = thermo.T()[celli];

        Info<< "t = " << t << ": <Y_A> " << moment("0_1_0")
            << " (exact " << decay*YA0[0] << "), <xi Y_A> "
            << moment("1_1_0") << ", T " << thermo.T()[celli] << " K";

        if (conditionalEnthalpy)
        {
            Info<< ", enthalpy off the line by " << worstLine
                << " in the nodes and " << worstMomentLine
                << " in the moments";
        }

        Info<< nl;
    }


    // * * * * * * * * * * * * * * * * The checks * * * * * * * * * * * * //

    Info<< nl << "Checks over " << nSteps << " steps:" << nl;

    check("the mean of Y_A decays at the rate", worstMeanDecay, 1e-8);

    if (!mixingOn)
    {
        check
        (
            "<xi Y_A> and <xi^2 Y_A> decay at the rate",
            worstDecay,
            1e-8
        );

        check
        (
            "the moments of the mixture fraction do not change",
            worstXi,
            1e-14
        );
    }

    check("the zero-order moment stays one", worstM0, 1e-12);
    check("the mean of the mixture fraction is conserved", worstMean, 1e-12);
    check
    (
        "<xi^k Y_A> + <xi^k Y_B> = <xi^k>, and every node sums to one",
        worstSum,
        1e-12
    );
    check("every mass fraction stays within [0, 1]", worstBound, 0);
    check("the thermo holds the mean composition", worstMeanY, 1e-12);
    check("every environment's field holds its temperature", worstNodeT, 1e-10);

    if (conditionalEnthalpy)
    {
        check
        (
            "the moments of the enthalpy stay linear in those of xi",
            worstMomentLine,
            1e-8
        );

        check
        (
            "the enthalpy of every environment stays linear in xi",
            worstLine,
            1e-5
        );
    }

    if (linearEnthalpy)
    {
        check
        (
            "the temperature of the cell is the mean of the environments",
            worstT,
            1e-8
        );

        check
        (
            "the temperature rises with the heat of reaction",
            Tinitial - Tfinal,
            -1
        );

        Info<< "    (from " << Tinitial << " to " << Tfinal << " K)" << nl;
    }
    else
    {
        check
        (
            "the temperature of the cell does not change",
            mag(Tfinal - Tinitial),
            1e-12
        );
    }

    Info<< nl;

    if (nFailed)
    {
        Info<< nFailed << " of " << nChecks << " checks FAILED" << nl << endl;
        return 1;
    }

    Info<< "All " << nChecks << " checks passed" << nl << endl;

    return 0;
}


// ************************************************************************* //
