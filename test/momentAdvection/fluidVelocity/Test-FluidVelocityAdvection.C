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
    Test-FluidVelocityAdvection

Description
    Test the transport of a quadrature of conditional means by the velocity
    of the fluid that carries it, on a periodic line where the answer is
    known: a uniform velocity moves the mixture fraction and the
    composition conditioned on it around the line without changing them.

    The distribution is built from three nodes of the mixture fraction, each
    carrying one composition, whose weights, abscissae and compositions vary
    from cell to cell. Its moments are set from those nodes, so that every
    check is against a known answer:

    - every moment is conserved, whatever the Courant number;
    - at a Courant number of one the moments come back exactly where they
      started after one period, since upwinding is a shift then;
    - the moment set of every cell stays realizable and every conditional
      mean stays within [0, 1], at the Courant number the scheme allows;
    - the answer does not know left from right: the mirror image of the
      profile, transported the other way, is the mirror image of the result.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "quadratureApproximations.H"
#include "velocityMomentAdvection.H"
#include "univariateMomentSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nChecks = 0;
label nFailed = 0;

void check(const word& name, const bool ok, const string& detail = "")
{
    nChecks++;

    if (ok)
    {
        Info<< "    passed: " << name << nl;
    }
    else
    {
        nFailed++;
        Info<< "    FAILED: " << name << " " << detail.c_str() << nl;
    }
}


//- The nodes of cell celli of nCells: three nodes of the mixture fraction,
//  each carrying a composition, in a profile that is not symmetric so that
//  the mirror check has something to catch
void nodesOf
(
    const label celli,
    const label nCells,
    const bool mirrored,
    scalarList& ws,
    scalarList& xs,
    scalarList& ys
)
{
    const label i = mirrored ? nCells - 1 - celli : celli;
    const scalar s = (scalar(i) + 0.5)/scalar(nCells);
    const scalar bump = Foam::exp(-30.0*sqr(s - 0.35));

    ws = scalarList({0.2 + 0.5*bump, 0.5, 0.3 + 0.2*s});
    xs = scalarList({0.1*(1 + 0.5*s), 0.45 + 0.2*bump, 0.9 - 0.05*s});
    ys = scalarList({0.05 + 0.9*bump, 0.5 - 0.3*s, 0.95});
}


//- Set the moments of the quadrature from those nodes
void setMoments
(
    volVelocityMomentFieldSet& moments,
    const labelListList& momentOrders,
    const bool mirrored
)
{
    const label nCells = moments(0).mesh().nCells();

    scalarList ws, xs, ys;

    for (label celli = 0; celli < nCells; celli++)
    {
        nodesOf(celli, nCells, mirrored, ws, xs, ys);

        forAll(momentOrders, mi)
        {
            const labelList& order = momentOrders[mi];
            scalar m = 0;

            forAll(ws, nodei)
            {
                m += ws[nodei]*pow(xs[nodei], order[0])
                    *pow(ys[nodei], order[1]);
            }

            moments[mi][celli] = m;
        }
    }

    forAll(momentOrders, mi)
    {
        moments[mi].correctBoundaryConditions();
    }
}


//- Whether the moments of the mixture fraction of every cell are
//  realizable on [0, 1], which is what the inversion tests
bool realizableEverywhere
(
    const volVelocityMomentFieldSet& moments,
    const labelListList& momentOrders
)
{
    const label nCells = moments(0).mesh().nCells();
    label nPrimary = 0;

    forAll(momentOrders, mi)
    {
        if (momentOrders[mi][1] == 0)
        {
            nPrimary++;
        }
    }

    for (label celli = 0; celli < nCells; celli++)
    {
        univariateMomentSet primary
        (
            nPrimary, supportType::ZeroOne, SMALL, 0, Zero
        );

        forAll(primary, mi)
        {
            primary[mi] = moments(labelList({mi, 0}))[celli];
        }

        if (!primary.isRealizable(false))
        {
            return false;
        }
    }

    return true;
}


//- The moments of the whole line, which transport conserves
scalarField totalMoments(const volVelocityMomentFieldSet& moments)
{
    scalarField total(moments.size(), Zero);

    forAll(total, mi)
    {
        total[mi] = gSum(moments[mi].primitiveField()*moments(0).mesh().V());
    }

    return total;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // The mixture fraction and the composition conditioned on it are both
    // on [0, 1]: the supports are the model's, and this test is the model
    const List<supportType> supports(2, supportType::ZeroOne);

    velocityQuadratureApproximation quadrature("mixture", mesh, supports);

    volVelocityMomentFieldSet& moments = quadrature.moments();
    const labelListList momentOrders(quadrature.momentOrders());
    const label nCells = mesh.nCells();
    const label nMoments = moments.size();

    const dictionary& schemeDict = quadrature.subDict("momentAdvection");

    // A uniform velocity along the line, of one cell per step at a Courant
    // number of one. The spacing comes from the mesh, and the flux of a
    // face from its area, so that the Courant number is the one meant
    // whatever the mesh is.
    const boundBox bounds(mesh.points());
    const scalar dx = (bounds.max().x() - bounds.min().x())/scalar(nCells);
    const scalar U0 = 1.0;

    Info<< "\nA quadrature of " << nMoments << " moments on "
        << nCells << " cells" << nl;

    //- Transport the profile, mirrored or not, for nSteps at a Courant
    //  number, checking what holds along the way, and return what every
    //  cell holds at the end
    auto advect =
        [&](const bool mirrored, const scalar Co, const label nSteps)
    {
        Info<< nl << (mirrored ? "Mirrored, " : "As written, ")
            << "Courant number " << Co << ", " << nSteps << " steps" << nl;

        const scalar deltaT = Co*dx/U0;

        surfaceScalarField phi
        (
            IOobject("phi", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("phi", dimVolume/dimTime, Zero)
        );

        const scalar U = mirrored ? -U0 : U0;
        const vector along(1, 0, 0);

        phi.primitiveFieldRef() = U*(mesh.Sf().primitiveField() & along);

        forAll(phi.boundaryField(), patchi)
        {
            phi.boundaryFieldRef()[patchi] =
                U*(mesh.Sf().boundaryField()[patchi] & along);
        }

        setMoments(moments, momentOrders, mirrored);
        quadrature.updateQuadrature();

        const scalarField initial(totalMoments(moments));

        autoPtr<velocityMomentAdvection> advection
        (
            velocityMomentAdvection::New(schemeDict, quadrature, supports)
        );

        check
        (
            "the scheme allows a Courant number of one",
            advection->realizableCo() >= 1.0
        );

        scalar worstConservation = 0;
        bool realizable = true;
        bool bounded = true;

        for (label stepi = 0; stepi < nSteps; stepi++)
        {
            advection->update();

            forAll(moments, mi)
            {
                moments[mi].primitiveFieldRef() -=
                    deltaT*advection->divMoments()[mi].primitiveField();

                moments[mi].correctBoundaryConditions();
            }

            realizable =
                realizable && realizableEverywhere(moments, momentOrders);

            quadrature.updateQuadrature();

            forAll(quadrature.nodes(), nodei)
            {
                const volVelocityNode& node = quadrature.nodes()[nodei];

                forAll(node.abscissae()[1], celli)
                {
                    const scalar y = node.abscissae()[1][celli];

                    bounded = bounded && y >= 0 && y <= 1;
                }
            }

            const scalarField total(totalMoments(moments));

            forAll(total, mi)
            {
                worstConservation =
                    max
                    (
                        worstConservation,
                        mag(total[mi] - initial[mi])
                       /max(mag(initial[mi]), SMALL)
                    );
            }
        }

        check
        (
            "every moment is conserved",
            worstConservation < 1e-12,
            "(" + Foam::name(worstConservation) + ")"
        );

        check("every cell stays realizable", realizable);
        check("every conditional mean stays within [0, 1]", bounded);

        scalarField result(nMoments*nCells);

        forAll(moments, mi)
        {
            forAll(moments[mi], celli)
            {
                result[mi*nCells + celli] = moments[mi][celli];
            }
        }

        return result;
    };

    // One period at a Courant number of one is a shift by the whole line,
    // which upwinding does exactly
    {
        setMoments(moments, momentOrders, false);

        scalarField start(nMoments*nCells);

        forAll(moments, mi)
        {
            forAll(moments[mi], celli)
            {
                start[mi*nCells + celli] = moments[mi][celli];
            }
        }

        const scalarField after(advect(false, 1.0, nCells));

        check
        (
            "one period at a Courant number of one is a shift",
            gMax(mag(after - start)) < 1e-12,
            "(" + Foam::name(gMax(mag(after - start))) + ")"
        );
    }

    // Half a line at half that step, where upwinding smears the profile:
    // the mirror image of the profile, carried the other way, has to give
    // the mirror image of the result
    {
        const scalarField forwards(advect(false, 0.5, nCells));
        const scalarField backwards(advect(true, 0.5, nCells));

        scalar worstMirror = 0;

        for (label mi = 0; mi < nMoments; mi++)
        {
            for (label celli = 0; celli < nCells; celli++)
            {
                worstMirror =
                    max
                    (
                        worstMirror,
                        mag
                        (
                            forwards[mi*nCells + celli]
                          - backwards[mi*nCells + nCells - 1 - celli]
                        )
                    );
            }
        }

        check
        (
            "the transport does not know left from right",
            worstMirror < 1e-12,
            "(" + Foam::name(worstMirror) + ")"
        );
    }

    Info<< nl;

    if (nFailed > 0)
    {
        FatalErrorInFunction
            << nFailed << " of " << nChecks << " checks failed." << nl
            << exit(FatalError);
    }

    Info<< "All " << nChecks << " checks passed" << nl << endl;

    return 0;
}

// ************************************************************************* //
