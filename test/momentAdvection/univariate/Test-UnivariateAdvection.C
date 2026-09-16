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
    Test-UnivariateAdvection

Description
    Checks the univariate moment advection schemes on a periodic
    one-dimensional mesh, where the transport of a moment set by a uniform
    velocity has properties that hold exactly.

    The scheme is driven directly, the way the PDF transport model drives
    it: the moments are updated by the divergence the scheme returns, with
    the Euler step the transport model takes and no source, so what is
    checked is the scheme alone. The profile is set by the test, from the
    same formulas as setMoments1D, and the velocity is uniform, so that
    after one period the population is back where it started.

    The checks are selected by the case, in system/advectionTestDict:

    - conservation: the integral of every moment over the domain does not
      change, step by step, whatever the Courant number.
    - realizability: every cell that holds a population holds a realizable
      moment set at every step. This is what the schemes guarantee, under
      the Courant limit each of them reports.
    - mirror: the profile advected to the right and its mirror image
      advected to the left are mirror images of one another at the end,
      to the tolerance the case sets: round-off for a scheme without a
      limiter, and for one whose limiter decides on a realizability test
      what the cells that test settles either way within round-off can
      move. A scheme that treats the two sides of a face differently
      fails it by orders of magnitude.
    - shift: at a Courant number of one the first-order scheme moves the
      profile by exactly one cell per step, so after a period the moments
      equal the initial ones to round-off.
    - accuracy: the scheme selected by the case, run over one period, is
      closer to the exact solution than the first-order scheme selected
      by momentAdvectionFirstOrder is, by the ratio the case sets.

    Run from the case directory.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "supportType.H"
#include "quadratureApproximations.H"
#include "univariateMomentAdvection.H"
#include "univariateMomentSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * The profiles  * * * * * * * * * * * * * * //

//- The moments of a beta distribution whose parameters vary along x,
//  ramped to zero at both ends of the domain: what setMoments1D calls
//  regular. Its support is [0, 1], so it serves both supports.
void regularProfile(const scalar x, scalarList& m)
{
    using constant::mathematical::pi;

    const scalar alpha = 3.5 + 1.5*Foam::sin(2.0*pi*x);
    const scalar beta = 3.5 - 1.5*Foam::cos(2.0*pi*x);

    m[0] = 1.0;

    for (label k = 1; k < m.size(); k++)
    {
        m[k] = m[k - 1]*(alpha + scalar(k - 1))/(alpha + beta + scalar(k - 1));
    }

    const scalar ramp =
        x <= 0.5
      ? 0.5*(1.0 + Foam::tanh(Foam::tan(pi*(2.0*x - 0.5))))
      : 0.5*(1.0 + Foam::tanh(Foam::tan(pi*(-2.0*x + 1.5))));

    forAll(m, k)
    {
        m[k] *= ramp;
    }
}


//- Two Dirac modes and a continuous one whose weights and parameters vary
//  along x, so that the number of nodes the moments hold changes from cell
//  to cell: what setMoments1D calls bimodal, and what makes the limiter of
//  the second-order scheme act. Its support is R+.
void bimodalProfile(const scalar x, scalarList& m)
{
    const scalar w1 = (x >= 0.0 && x <= 1.0) ? 16.0*sqr(x)*sqr(1.0 - x) : 0.0;

    const scalar w2 =
        (4.0*x >= 1.0 && x <= 1.0)
      ? 256.0*sqr(4.0*x - 1.0)*sqr(1.0 - x)/81.0
      : 0.0;

    const scalar w3 =
        (3.0*x >= 1.0 && x <= 1.0) ? 9.0*sqr(3.0*x - 1.0)*sqr(1.0 - x) : 0.0;

    auto blend = [x](const scalar lower, const scalar upper)
    {
        if (3.0*x <= 1.0)
        {
            return lower;
        }
        else if (3.0*x <= 2.0)
        {
            return
                lower*sqr(2.0 - 3.0*x)*(6.0*x - 1.0)
              + upper*sqr(3.0*x - 1.0)*(5.0 - 6.0*x);
        }

        return upper;
    };

    const scalar lambda2 = blend(2.0e-2, 0.7);
    const scalar k2 = blend(3.0, 10.0);
    const scalar xi1 = 0.02;

    forAll(m, k)
    {
        m[k] =
            w1*Foam::pow(xi1, k)
          + w2*Foam::pow(2.0*xi1, k)
          + w3*Foam::pow(lambda2, k)*std::tgamma(1.0 + scalar(k)/k2);
    }
}


// * * * * * * * * * * * * * * * * The driver * * * * * * * * * * * * * * * //

//- What a run needs to know about the mesh: the cells sorted along x, so
//  that a cell and its mirror image can be paired, and the spacing
struct meshLine
{
    labelList sorted;
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
    }

    //- The cell that mirrors a cell about the middle of the domain
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


//- Set the moments of every cell from the profile, as it is or mirrored
//  about the middle of the domain
void setProfile
(
    const word& profile,
    const bool mirrored,
    const fvMesh& mesh,
    const meshLine& mesh1D,
    volScalarMomentFieldSet& moments
)
{
    scalarList m(moments.size(), Zero);

    forAll(mesh.C(), celli)
    {
        // The domain is the unit interval: the position is measured from
        // its lower end, whatever the mesh coordinates are
        scalar x =
            (mesh.C()[celli].x() - mesh.C()[mesh1D.sorted[0]].x())
           /mesh1D.length
          + 0.5/scalar(mesh.nCells());

        if (mirrored)
        {
            x = 1.0 - x;
        }

        if (profile == "regular")
        {
            regularProfile(x, m);
        }
        else if (profile == "bimodal")
        {
            bimodalProfile(x, m);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown profile " << profile
                << ": valid profiles are regular and bimodal"
                << exit(FatalError);
        }

        forAll(moments, mi)
        {
            moments[mi][celli] = m[mi];
        }
    }

    forAll(moments, mi)
    {
        moments[mi].correctBoundaryConditions();
    }
}


//- The integral of every moment over the domain
scalarList integrals(const fvMesh& mesh, const volScalarMomentFieldSet& moments)
{
    scalarList sums(moments.size(), Zero);

    forAll(moments, mi)
    {
        sums[mi] = gSum(moments[mi].primitiveField()*mesh.V());
    }

    return sums;
}


//- The largest magnitude of every moment over the domain, the scale its
//  differences are measured against
scalarList scales(const volScalarMomentFieldSet& moments)
{
    scalarList s(moments.size(), Zero);

    forAll(moments, mi)
    {
        s[mi] = gMax(mag(moments[mi].primitiveField()));
    }

    return s;
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

    const word profile(testDict.get<word>("profile"));
    const supportType support(wordToSupportType(testDict.get<word>("support")));
    const scalar Co = testDict.get<scalar>("Co");
    const wordList checks(testDict.get<wordList>("checks"));

    auto wanted = [&checks](const word& name)
    {
        return checks.found(name);
    };

    scalarQuadratureApproximation quadrature
    (
        "populationBalance",
        mesh,
        List<supportType>(1, support)
    );

    volScalarMomentFieldSet& moments = quadrature.moments();
    const label nMoments = moments.size();

    const meshLine mesh1D(mesh);
    const label nCells = mesh.nCells();

    // One period of the uniform velocity brings the population back to
    // where it started; the number of steps in it is the number of cells
    // over the Courant number, which the case has to make an integer
    const scalar U0 = 1.0;
    scalar deltaT = Co*mesh1D.dx/U0;
    const label nSteps = label(scalar(nCells)/Co + 0.5);

    Info<< "\nAdvecting the " << profile << " profile over " << nCells
        << " cells at Co = " << Co << ", " << nSteps << " steps per period"
        << nl << endl;

    // The tolerances of the realizability check: a cell holding less than
    // a vanishing part of the largest population holds none, and a zeta
    // below round-off is on the boundary of the moment space
    const scalar smallZeta = 1.0e-10;

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


    // * * * * * * * * * * * * * * * A run * * * * * * * * * * * * * * * * //

    //- Advect the profile for a number of steps with the scheme of a
    //  dictionary, checking conservation and realizability as it goes when
    //  asked to. The velocity is along +x for a positive sign.
    auto run = [&]
    (
        const dictionary& schemeDict,
        const scalar direction,
        const label steps,
        const bool checkAsItGoes,
        const word& runName
    )
    {
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
            mesh.Sf() & dimensionedVector("U", dimVelocity, vector(direction*U0, 0, 0))
        );

        autoPtr<univariateMomentAdvection> advection
        (
            univariateMomentAdvection::New(schemeDict, quadrature, phi, support)
        );

        runTime.setDeltaT(deltaT);

        const scalarList initial(integrals(mesh, moments));
        const scalarList scale(scales(moments));
        const scalar smallM0 = 1.0e-12*scale[0];

        scalar worstConservation = 0;
        bool realizable = true;
        label firstUnrealizableStep = -1;
        label firstUnrealizableCell = -1;

        for (label stepi = 0; stepi < steps; stepi++)
        {
            runTime++;

            advection().update();

            const mappedPtrList<volScalarField>& div = advection().divMoments();

            forAll(moments, mi)
            {
                moments[mi].primitiveFieldRef() -=
                    deltaT*div[mi].primitiveField();

                moments[mi].correctBoundaryConditions();
            }

            if (!checkAsItGoes)
            {
                continue;
            }

            const scalarList now(integrals(mesh, moments));

            forAll(now, mi)
            {
                worstConservation =
                    max
                    (
                        worstConservation,
                        mag(now[mi] - initial[mi])/max(mag(initial[mi]), VSMALL)
                    );
            }

            if (realizable)
            {
                scalarList m(nMoments, Zero);

                forAll(mesh.C(), celli)
                {
                    if (moments[0][celli] < smallM0)
                    {
                        continue;
                    }

                    forAll(m, mi)
                    {
                        m[mi] = moments[mi][celli];
                    }

                    univariateMomentSet cellMoments(m, support, smallM0, smallZeta);

                    if (!cellMoments.isRealizable(false))
                    {
                        realizable = false;
                        firstUnrealizableStep = stepi;
                        firstUnrealizableCell = celli;
                        break;
                    }
                }
            }
        }

        if (checkAsItGoes)
        {
            check
            (
                worstConservation < 1.0e-12,
                runName + ": every moment is conserved to round-off (worst "
              + Foam::name(worstConservation) + ")"
            );

            check
            (
                realizable,
                runName + ": every populated cell stays realizable"
              + (
                    realizable
                  ? string("")
                  : string(" (first failure at step ")
                  + Foam::name(firstUnrealizableStep) + ", cell "
                  + Foam::name(firstUnrealizableCell) + ")"
                )
            );
        }
    };


    // * * * * * * * * * * * * * * * The checks * * * * * * * * * * * * * //

    const dictionary& schemeDict = quadrature.subDict("momentAdvection");

    // The scheme's own Courant limit has to admit the step the case asks
    // for, or the realizability it guarantees is not guaranteed
    setProfile(profile, false, mesh, mesh1D, moments);

    {
        surfaceScalarField phi
        (
            IOobject("phi", runTime.timeName(), mesh),
            mesh.Sf() & dimensionedVector("U", dimVelocity, vector(U0, 0, 0))
        );

        autoPtr<univariateMomentAdvection> advection
        (
            univariateMomentAdvection::New(schemeDict, quadrature, phi, support)
        );

        Info<< "The scheme allows Co up to " << advection().realizableCo()
            << nl << endl;

        check
        (
            Co <= advection().realizableCo() + SMALL,
            "the Courant number of the case is within the scheme's limit"
        );
    }

    // A forward run over one period, with conservation and realizability
    // checked at every step
    setProfile(profile, false, mesh, mesh1D, moments);

    PtrList<scalarField> initialMoments(nMoments);

    forAll(moments, mi)
    {
        initialMoments.set(mi, new scalarField(moments[mi].primitiveField()));
    }

    run
    (
        schemeDict,
        1.0,
        nSteps,
        wanted("conservation") || wanted("realizability"),
        "forward"
    );

    PtrList<scalarField> forwardMoments(nMoments);

    forAll(moments, mi)
    {
        forwardMoments.set(mi, new scalarField(moments[mi].primitiveField()));
    }

    const scalarList scale(scales(moments));

    if (wanted("shift"))
    {
        // At Co = 1 the first-order scheme moves every cell's moments to
        // the next cell each step, so a period brings them back exactly
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
            "a period at Co = 1 returns the moments to round-off (worst "
          + Foam::name(worst) + ")"
        );
    }

    if (wanted("mirror"))
    {
        // The mirrored profile advected to the left has to be the mirror
        // image of the profile advected to the right
        setProfile(profile, true, mesh, mesh1D, moments);

        run
        (
            schemeDict,
            -1.0,
            nSteps,
            wanted("conservation") || wanted("realizability"),
            "backward"
        );

        scalar worst = 0;
        label worstCell = -1;

        forAll(moments, mi)
        {
            forAll(mesh.C(), celli)
            {
                const scalar difference =
                    mag
                    (
                        moments[mi][mesh1D.mirror(celli)]
                      - forwardMoments[mi][celli]
                    )/scale[mi];

                if (difference > worst)
                {
                    worst = difference;
                    worstCell = celli;
                }
            }
        }

        // Where the two differ, for the record
        {
            DynamicList<label> differing;

            forAll(mesh.C(), celli)
            {
                forAll(moments, mi)
                {
                    const scalar difference =
                        mag
                        (
                            moments[mi][mesh1D.mirror(celli)]
                          - forwardMoments[mi][celli]
                        )/scale[mi];

                    if (difference > 1.0e-10)
                    {
                        differing.append(celli);
                        break;
                    }
                }
            }

            Info<< "  mirror: " << differing.size()
                << " cells differ by more than 1e-10" << endl;
        }

        // The tolerance is the case's: a scheme without a limiter is
        // symmetric to round-off, while the second-order scheme decides
        // whether to limit a cell on a realizability test of m*, which the
        // two runs can settle differently where m* lies on the boundary of
        // the moment space within round-off. Those cells cost little, so
        // the runs still agree to a tolerance a scheme that treats the two
        // sides of a face differently misses by orders of magnitude.
        const scalar tolerance = testDict.get<scalar>("mirrorTolerance");

        check
        (
            worst < tolerance,
            "the profile advected left is the mirror image of the profile "
            "advected right to " + Foam::name(tolerance) + " (worst "
          + Foam::name(worst) + " at cell " + Foam::name(worstCell) + ")"
        );
    }

    if (wanted("accuracy"))
    {
        // After a period the exact solution is the initial profile: the
        // scheme has to beat the first-order one by the ratio the case
        // sets. This runs at its own Courant number: the transport model
        // takes an Euler step, whose error grows with the step, so at the
        // Courant limit of the scheme the time error hides the order of
        // the reconstruction, and only a small step lets it show.
        const scalar ratio = testDict.get<scalar>("accuracyRatio");
        const scalar accuracyCo = testDict.get<scalar>("accuracyCo");

        deltaT = accuracyCo*mesh1D.dx/U0;
        const label accuracySteps = label(scalar(nCells)/accuracyCo + 0.5);

        Info<< nl << "Comparing with the first-order scheme at Co = "
            << accuracyCo << ", " << accuracySteps << " steps per period"
            << endl;

        setProfile(profile, false, mesh, mesh1D, moments);
        run(schemeDict, 1.0, accuracySteps, false, "scheme");

        PtrList<scalarField> schemeMoments(nMoments);

        forAll(moments, mi)
        {
            schemeMoments.set(mi, new scalarField(moments[mi].primitiveField()));
        }

        setProfile(profile, false, mesh, mesh1D, moments);

        run
        (
            quadrature.subDict("momentAdvectionFirstOrder"),
            1.0,
            accuracySteps,
            false,
            "first order"
        );

        forAll(moments, mi)
        {
            const scalar errorScheme =
                gSum(mag(schemeMoments[mi] - initialMoments[mi]));

            const scalar errorFirstOrder =
                gSum(mag(moments[mi].primitiveField() - initialMoments[mi]));

            check
            (
                errorScheme < ratio*errorFirstOrder,
                "moment " + Foam::name(mi) + ": the scheme's error is below "
              + Foam::name(ratio) + " of the first-order one ("
              + Foam::name(errorScheme/errorFirstOrder) + ")"
            );
        }
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
