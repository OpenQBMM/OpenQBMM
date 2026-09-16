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
    Test-VelocityAdvection

Description
    Checks the velocity moment advection schemes and their boundary
    conditions on a one-dimensional mesh, with populations whose transport
    has properties that hold exactly.

    The scheme is driven directly, the way the PDF transport model drives
    it: the quadrature is inverted, the scheme updated, and the moments
    moved by the divergence it returns with the Euler step the transport
    model takes and no source, so what is checked is the scheme and its
    boundary conditions alone.

    The population is a sum of beams the case describes, each a Gaussian
    in velocity with a mean, a temperature in each direction (zero makes
    the beam monokinetic) and a volume fraction that is a smooth bump
    along x, so that every moment is known in closed form at the start and
    after any translation.

    The checks are selected by the case, in system/advectionTestDict:

    - conservation: the volume fraction, the momentum and the kinetic
      energy of the population do not change on a periodic mesh, whatever
      the Courant number; with walls, the volume fraction only.
    - realizability: every populated cell keeps a non-negative volume
      fraction and non-negative variances in every direction.
    - shift: at a Courant number of one, monokinetic beams of equal speed
      are moved by exactly one cell per step, so after a period the
      moments equal the initial ones to round-off — beams crossing at
      opposite velocities pass through one another exactly.
    - mirror: the beams advected to the right and their mirror images
      advected to the left are mirror images of one another at the end.
    - energy: with walls, the kinetic energy of the population does not
      grow; with a restitution coefficient of one it is conserved.
    - bulkVelocity: with walls, the mean velocity along the wall of the
      population does not change.
    - heating: with walls hotter than the population, its temperature
      rises.
    - emptied: with an outflow, a beam that has had the time to leave has
      left, and nothing came back.

    Run from the case directory.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "supportType.H"
#include "quadratureApproximations.H"
#include "velocityMomentAdvection.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * The beams * * * * * * * * * * * * * * * //

//- A beam: a Gaussian in velocity, spread along x as a smooth bump
struct beam
{
    scalar centre;
    scalar width;
    scalar alpha;
    vector U;
    vector Theta;

    beam()
    :
        centre(0), width(0), alpha(0), U(Zero), Theta(Zero)
    {}

    beam(const dictionary& dict)
    :
        centre(dict.get<scalar>("centre")),
        width(dict.get<scalar>("width")),
        alpha(dict.get<scalar>("alpha")),
        U(dict.get<vector>("U")),
        Theta(dict.get<vector>("Theta"))
    {}

    //- The volume fraction at x on the unit periodic interval: a cosine
    //  bump of the given width about the centre, zero outside it
    scalar fraction(const scalar x) const
    {
        scalar d = x - centre;
        d -= std::round(d);

        if (mag(d) >= 0.5*width)
        {
            return 0;
        }

        return alpha*sqr(Foam::cos(constant::mathematical::pi*d/width));
    }

    //- The raw moment of order n of a Gaussian of mean m and variance s
    static scalar gaussianMoment(const label n, const scalar m, const scalar s)
    {
        switch (n)
        {
            case 0: return 1;
            case 1: return m;
            case 2: return sqr(m) + s;
            case 3: return pow3(m) + 3*m*s;
            case 4: return pow4(m) + 6*sqr(m)*s + 3*sqr(s);
        }

        FatalErrorInFunction
            << "No Gaussian moment of order " << n << exit(FatalError);

        return 0;
    }

    //- The moment of the beam at x, of the given orders in u, v and w
    scalar moment(const scalar x, const labelList& order) const
    {
        scalar m = fraction(x);

        forAll(order, cmpt)
        {
            m *= gaussianMoment(order[cmpt], U[cmpt], Theta[cmpt]);
        }

        return m;
    }

    //- The mirror image of the beam about the middle of the interval
    beam mirrored() const
    {
        beam b(*this);
        b.centre = 1.0 - centre;
        b.U.x() = -U.x();
        return b;
    }
};


// * * * * * * * * * * * * * * * * The driver * * * * * * * * * * * * * * * //

//- The cells sorted along x, the spacing and the length of the line
struct meshLine
{
    labelList sorted;
    scalar x0;
    scalar dx;
    scalar length;

    meshLine(const fvMesh& mesh)
    {
        const label n = mesh.nCells();
        sorted = identity(n);

        const vectorField& C = mesh.C();

        std::sort
        (
            sorted.begin(),
            sorted.end(),
            [&C](label a, label b){ return C[a].x() < C[b].x(); }
        );

        dx = C[sorted[1]].x() - C[sorted[0]].x();
        length = scalar(n)*dx;
        x0 = C[sorted[0]].x() - 0.5*dx;
    }

    //- The position of a cell on the unit interval
    scalar x(const fvMesh& mesh, const label celli) const
    {
        return (mesh.C()[celli].x() - x0)/length;
    }

    //- The cell that mirrors a cell about the middle of the line
    label mirror(const label celli) const
    {
        const label n = sorted.size();

        forAll(sorted, i)
        {
            if (sorted[i] == celli)
            {
                return sorted[n - 1 - i];
            }
        }

        return -1;
    }
};


//- Set every moment of every cell from the beams
void setBeams
(
    const List<beam>& beams,
    const fvMesh& mesh,
    const meshLine& mesh1D,
    const labelListList& orders,
    volVelocityMomentFieldSet& moments
)
{
    forAll(mesh.C(), celli)
    {
        const scalar x = mesh1D.x(mesh, celli);

        forAll(moments, mi)
        {
            scalar m = 0;

            forAll(beams, bi)
            {
                m += beams[bi].moment(x, orders[mi]);
            }

            moments[mi][celli] = m;
        }
    }

    forAll(moments, mi)
    {
        moments[mi].correctBoundaryConditions();
    }
}


//- The index of a moment order in the set, or -1
label indexOf(const labelListList& orders, const labelList& order)
{
    forAll(orders, mi)
    {
        if (orders[mi] == order)
        {
            return mi;
        }
    }

    return -1;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary testDict
    (
        IOobject
        (
            "advectionTestDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar Co = testDict.get<scalar>("Co");
    const wordList checks(testDict.get<wordList>("checks"));

    auto wanted = [&checks](const word& name)
    {
        return checks.found(name);
    };

    List<beam> beams;

    {
        const List<dictionary> beamDicts(testDict.lookup("beams"));

        forAll(beamDicts, bi)
        {
            beams.append(beam(beamDicts[bi]));
        }
    }

    velocityQuadratureApproximation quadrature
    (
        "particles",
        mesh,
        List<supportType>(3, supportType::R)
    );

    volVelocityMomentFieldSet& moments = quadrature.moments();
    const labelListList& orders = quadrature.momentOrders();
    const label nMoments = moments.size();

    // The moments the checks read
    const label m000 = indexOf(orders, {0, 0, 0});
    const label m100 = indexOf(orders, {1, 0, 0});
    const label m010 = indexOf(orders, {0, 1, 0});
    const label m001 = indexOf(orders, {0, 0, 1});
    const label m200 = indexOf(orders, {2, 0, 0});
    const label m020 = indexOf(orders, {0, 2, 0});
    const label m002 = indexOf(orders, {0, 0, 2});

    const meshLine mesh1D(mesh);
    const label nCells = mesh.nCells();

    // The velocity scale the checks measure against: the fastest beam,
    // with the spread its temperature gives it
    scalar Umax = 0;

    forAll(beams, bi)
    {
        Umax =
            max
            (
                Umax,
                mag(beams[bi].U) + Foam::sqrt(3.0*cmptMax(beams[bi].Theta))
            );
    }

    // A period is the number of cells over the Courant number: for
    // monokinetic beams of equal speed that brings them back where they
    // started, and for any other population it is a length of time
    const scalar periods = testDict.lookupOrDefault<scalar>("periods", 1);
    const label nSteps = label(periods*scalar(nCells)/Co + 0.5);

    // The step, set below from the Courant number of the scheme
    scalar deltaT = 1;

    Info<< "\nAdvecting " << beams.size() << " beam(s) over " << nCells
        << " cells at Co = " << Co << ", " << nSteps << " steps per period"
        << nl << endl;

    label nTested = 0;
    label nFailed = 0;

    auto check = [&](const bool passed, const string& what)
    {
        nTested++;

        if (passed)
        {
            Info<< "  passed: " << what.c_str() << endl;
        }
        else
        {
            nFailed++;
            Info<< "  FAILED: " << what.c_str() << endl;
        }
    };

    const dictionary& schemeDict = quadrature.subDict("momentAdvection");

    // The integrals the conservation checks watch: volume fraction, the
    // three components of momentum, and twice the kinetic energy
    auto integrals = [&]()
    {
        scalarList sums(5, Zero);

        forAll(mesh.V(), celli)
        {
            const scalar V = mesh.V()[celli];

            sums[0] += moments[m000][celli]*V;
            sums[1] += moments[m100][celli]*V;
            sums[2] += moments[m010][celli]*V;
            sums[3] += moments[m001][celli]*V;
            sums[4] +=
                (
                    moments[m200][celli]
                  + moments[m020][celli]
                  + moments[m002][celli]
                )*V;
        }

        return sums;
    };

    //- Advect the beams for a number of steps, watching the integrals and
    //  the realizability as it goes. Returns the largest relative drift
    //  of each integral, and the growth of the energy.
    auto run = [&]
    (
        const label steps,
        const word& runName,
        scalarList& drift,
        scalar& energyGrowth,
        bool& realizable
    )
    {
        autoPtr<velocityMomentAdvection> advection
        (
            velocityMomentAdvection::New
            (
                schemeDict,
                quadrature,
                List<supportType>(3, supportType::R)
            )
        );

        runTime.setDeltaT(deltaT);

        const scalarList initial(integrals());
        scalar alphaScale = 0;

        forAll(mesh.C(), celli)
        {
            alphaScale = max(alphaScale, moments[m000][celli]);
        }

        drift = scalarList(5, Zero);
        energyGrowth = 0;
        realizable = true;

        for (label stepi = 0; stepi < steps; stepi++)
        {
            runTime++;

            quadrature.updateQuadrature();
            advection().update();

            const PtrList<volScalarField>& div = advection().divMoments();

            forAll(moments, mi)
            {
                moments[mi].primitiveFieldRef() -=
                    deltaT*div[mi].primitiveField();

                moments[mi].correctBoundaryConditions();
            }

            const scalarList now(integrals());

            forAll(now, i)
            {
                // The momentum along a direction the beams do not move in
                // is zero, and its drift is measured against the volume
                // fraction times the fastest speed instead
                const scalar scale =
                    i == 0 ? mag(initial[0])
                  : i < 4 ? max(mag(initial[i]), initial[0]*Umax)
                  : mag(initial[4]);

                drift[i] = max(drift[i], mag(now[i] - initial[i])/scale);
            }

            energyGrowth = max(energyGrowth, (now[4] - initial[4])/initial[4]);

            if (realizable)
            {
                forAll(mesh.C(), celli)
                {
                    const scalar a = moments[m000][celli];

                    if (a < -1.0e-12*alphaScale)
                    {
                        realizable = false;
                        break;
                    }

                    if (a < 1.0e-8*alphaScale)
                    {
                        continue;
                    }

                    const scalar Ux = moments[m100][celli]/a;
                    const scalar Uy = moments[m010][celli]/a;
                    const scalar Uz = moments[m001][celli]/a;

                    const scalar varianceScale = sqr(Umax);

                    if
                    (
                        moments[m200][celli]/a - sqr(Ux) < -1.0e-10*varianceScale
                     || moments[m020][celli]/a - sqr(Uy) < -1.0e-10*varianceScale
                     || moments[m002][celli]/a - sqr(Uz) < -1.0e-10*varianceScale
                    )
                    {
                        realizable = false;
                        break;
                    }
                }
            }
        }

        Info<< "  " << runName << ": drift of volume fraction "
            << drift[0] << ", momentum (" << drift[1] << " " << drift[2]
            << " " << drift[3] << "), energy " << drift[4] << endl;
    };


    // * * * * * * * * * * * * * * * The checks * * * * * * * * * * * * * //

    setBeams(beams, mesh, mesh1D, orders, moments);

    // The step is set from the Courant number the scheme itself reports
    // for a unit step, which is built on the velocity abscissae of the
    // nodes rather than on the mean velocity of a beam: a thermal beam
    // carries nodes well faster than its mean, and the kinetic scheme is
    // stable only while none of them crosses more than a cell per step.
    // For monokinetic beams the two coincide, and a period is a shift by
    // the number of cells.
    quadrature.updateQuadrature();

    {
        autoPtr<velocityMomentAdvection> advection
        (
            velocityMomentAdvection::New
            (
                schemeDict,
                quadrature,
                List<supportType>(3, supportType::R)
            )
        );

        runTime.setDeltaT(1.0);
        deltaT = Co/advection().CoNum();
    }

    Info<< "Step " << deltaT << " for Co = " << Co
        << " on the fastest node" << nl << endl;

    // The step has to be within the Courant limit the scheme reports for
    // it, or what the scheme guarantees is not guaranteed
    {
        autoPtr<velocityMomentAdvection> advection
        (
            velocityMomentAdvection::New
            (
                schemeDict,
                quadrature,
                List<supportType>(3, supportType::R)
            )
        );

        runTime.setDeltaT(deltaT);

        const scalar limit = advection().realizableCo();

        Info<< "The scheme allows a step up to " << limit
            << " times this one" << nl << endl;

        check
        (
            limit >= 1.0 - SMALL,
            "the step of the case is within the scheme's Courant limit ("
          + Foam::name(limit) + " times the step)"
        );
    }

    setBeams(beams, mesh, mesh1D, orders, moments);

    PtrList<scalarField> initialMoments(nMoments);

    forAll(moments, mi)
    {
        initialMoments.set(mi, new scalarField(moments[mi].primitiveField()));
    }

    scalarList drift;
    scalar energyGrowth;
    bool realizable;

    run(nSteps, "forward", drift, energyGrowth, realizable);

    PtrList<scalarField> forwardMoments(nMoments);

    forAll(moments, mi)
    {
        forwardMoments.set(mi, new scalarField(moments[mi].primitiveField()));
    }

    // The scale each moment is measured against: what it is, or what a
    // moment of its order would be for the population, whichever is
    // larger. A moment odd in a direction the population is symmetric in
    // starts at zero and carries round-off, and is measured against the
    // latter.
    scalarList scale(nMoments, Zero);

    {
        scalar alphaMax = gMax(initialMoments[m000]);

        forAll(moments, mi)
        {
            label order = 0;

            forAll(orders[mi], cmpt)
            {
                order += orders[mi][cmpt];
            }

            scale[mi] =
                max
                (
                    gMax(mag(initialMoments[mi])),
                    alphaMax*Foam::pow(Umax, order)
                );
        }
    }

    const scalar conservationTolerance =
        testDict.lookupOrDefault<scalar>("conservationTolerance", 1.0e-12);

    if (wanted("conservation"))
    {
        check
        (
            drift[0] < conservationTolerance,
            "the volume fraction is conserved to round-off (drift "
          + Foam::name(drift[0]) + ")"
        );

        check
        (
            max(drift[1], max(drift[2], drift[3])) < conservationTolerance,
            "the momentum is conserved to round-off (drift "
          + Foam::name(max(drift[1], max(drift[2], drift[3]))) + ")"
        );

        check
        (
            drift[4] < conservationTolerance,
            "the kinetic energy is conserved to round-off (drift "
          + Foam::name(drift[4]) + ")"
        );
    }

    if (wanted("mass"))
    {
        check
        (
            drift[0] < 1.0e-12,
            "the volume fraction is conserved to round-off (drift "
          + Foam::name(drift[0]) + ")"
        );
    }

    if (wanted("energy"))
    {
        const scalar tolerance =
            testDict.lookupOrDefault<scalar>("energyTolerance", 1.0e-12);

        check
        (
            energyGrowth < tolerance,
            "the kinetic energy does not grow (largest growth "
          + Foam::name(energyGrowth) + ")"
        );

        if (testDict.lookupOrDefault<bool>("energyConserved", false))
        {
            check
            (
                drift[4] < 1.0e-12,
                "the kinetic energy is conserved to round-off (drift "
              + Foam::name(drift[4]) + ")"
            );
        }
    }

    if (wanted("bulkVelocity"))
    {
        // The mean velocity along the walls, over the population
        const scalarList initial
        {
            gSum(initialMoments[m010])/gSum(initialMoments[m000]),
            gSum(initialMoments[m001])/gSum(initialMoments[m000])
        };

        const scalarList final
        {
            gSum(forwardMoments[m010])/gSum(forwardMoments[m000]),
            gSum(forwardMoments[m001])/gSum(forwardMoments[m000])
        };

        const scalar change =
            max(mag(final[0] - initial[0]), mag(final[1] - initial[1]))/Umax;

        check
        (
            change < 1.0e-12,
            "the bulk velocity along the walls is conserved (change "
          + Foam::name(change) + " of the beam speed)"
        );
    }

    if (wanted("heating"))
    {
        auto temperature = [&](const PtrList<scalarField>& m)
        {
            const scalar a = gSum(m[m000]);

            return
                (gSum(m[m200]) + gSum(m[m020]) + gSum(m[m002]))/a
              - (
                    sqr(gSum(m[m100])) + sqr(gSum(m[m010]))
                  + sqr(gSum(m[m001]))
                )/sqr(a);
        };

        const scalar before = temperature(initialMoments);
        const scalar after = temperature(forwardMoments);

        check
        (
            after > before,
            "the walls heat the population (temperature "
          + Foam::name(before) + " to " + Foam::name(after) + ")"
        );
    }

    if (wanted("emptied"))
    {
        const scalar left = gSum(forwardMoments[m000])/gSum(initialMoments[m000]);
        const scalar lowest = gMin(forwardMoments[m000])/gMax(initialMoments[m000]);

        // Where what is left sits, for the record
        label fullest = 0;

        forAll(forwardMoments[m000], celli)
        {
            if (forwardMoments[m000][celli] > forwardMoments[m000][fullest])
            {
                fullest = celli;
            }
        }

        Info<< "  most of what is left is in cell " << fullest << " at x = "
            << mesh1D.x(mesh, fullest) << ": "
            << forwardMoments[m000][fullest]/gMax(initialMoments[m000])
            << " of the initial peak" << endl;

        // What the scheme leaves behind: the tail of an upwind scheme
        // never quite reaches zero
        const scalar tolerance =
            testDict.lookupOrDefault<scalar>("emptiedTolerance", 1.0e-12);

        check
        (
            left < tolerance && lowest > -1.0e-12,
            "the beam has left through the outflow and nothing came back "
            "(fraction left " + Foam::name(left) + ", allowed "
          + Foam::name(tolerance) + ")"
        );
    }

    if (wanted("realizability"))
    {
        check
        (
            realizable,
            "every populated cell keeps non-negative volume fraction and "
            "variances"
        );
    }

    if (wanted("shift"))
    {
        scalar worst = 0;

        forAll(moments, mi)
        {
            worst =
                max
                (
                    worst,
                    gMax(mag(forwardMoments[mi] - initialMoments[mi]))/scale[mi]
                );
        }

        check
        (
            worst < 1.0e-12,
            "a period at Co = 1 returns every moment to round-off (worst "
          + Foam::name(worst) + ")"
        );
    }

    if (wanted("mirror"))
    {
        List<beam> mirrored(beams);

        forAll(mirrored, bi)
        {
            mirrored[bi] = beams[bi].mirrored();
        }

        setBeams(mirrored, mesh, mesh1D, orders, moments);

        run(nSteps, "backward", drift, energyGrowth, realizable);

        // A moment odd in u changes sign under the mirror
        scalar worst = 0;

        forAll(moments, mi)
        {
            const scalar sign = orders[mi][0] % 2 == 0 ? 1.0 : -1.0;

            forAll(mesh.C(), celli)
            {
                worst =
                    max
                    (
                        worst,
                        mag
                        (
                            sign*moments[mi][mesh1D.mirror(celli)]
                          - forwardMoments[mi][celli]
                        )/scale[mi]
                    );
            }
        }

        const scalar tolerance =
            testDict.lookupOrDefault<scalar>("mirrorTolerance", 1.0e-12);

        check
        (
            worst < tolerance,
            "the beams advected left are the mirror image of the beams "
            "advected right to " + Foam::name(tolerance) + " (worst "
          + Foam::name(worst) + ")"
        );
    }

    Info<< nl;

    if (nFailed > 0)
    {
        FatalErrorInFunction
            << nFailed << " of " << nTested << " checks failed." << nl
            << exit(FatalError);
    }

    Info<< nTested << " checks passed." << nl << endl;

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
