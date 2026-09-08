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
    Test-MomentAdvection

Description
    Test the moment reconstruction of the zeta advection scheme.

    The scheme does not transport moments directly. It maps a moment set onto
    the zero-order moment and a chain of auxiliary quantities - the zeta chain
    for a measure with support over R+, the canonical moments for one with
    support over [0, 1] - reconstructs those on the faces, and maps them back.
    The map back is zetaToMoments and canonicalMomentsToMoments, and this is
    what is tested here: whatever else the reconstruction does, the map has to
    return the moment set it was given when it is handed that set's own
    auxiliary quantities.

    That identity is what makes the scheme conservative, and it is exercised
    once per face per timestep by the two tutorials that select zeta, so an
    error in it is an error in every result the scheme produces.

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "scalarList.H"
#include "supportType.H"
#include "univariateMomentSet.H"
#include "zetaUnivariateAdvection.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label nTested = 0;

void check(const bool condition, const string& what)
{
    nTested++;

    if (!condition)
    {
        FatalErrorInFunction
            << "Failed: " << what << nl
            << exit(FatalError);
    }

    Info<< "  OK: " << what << endl;
}


//- Compare two moment sets to a relative tolerance, scaled by the magnitude
//  of the moment so that the high-order moments, which are large, are not
//  held to an absolute tolerance they cannot meet
void checkMoments
(
    const scalarList& computed,
    const scalarList& expected,
    const scalar tolerance,
    const string& what
)
{
    nTested++;

    scalar worst = 0;
    label worstOrder = -1;

    forAll(expected, mi)
    {
        const scalar error =
            mag(computed[mi] - expected[mi])/max(mag(expected[mi]), SMALL);

        if (error > worst)
        {
            worst = error;
            worstOrder = mi;
        }
    }

    if (worst > tolerance)
    {
        FatalErrorInFunction
            << "Failed: " << what << nl
            << "    Largest relative error " << worst
            << " on the moment of order " << worstOrder << nl
            << "    Computed " << computed << nl
            << "    Expected " << expected << nl
            << exit(FatalError);
    }

    Info<< "  OK: " << what
        << " (largest relative error " << worst << ")" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Build a moment set from a quadrature, which is realizable by construction,
//  and compute the auxiliary quantities of it. The realizability check has to
//  be run before the zeta chain is read, because it is what fills it.
autoPtr<univariateMomentSet> momentSetFromQuadrature
(
    const scalarList& weights,
    const scalarList& abscissae,
    const label nMoments,
    const supportType& support
)
{
    autoPtr<univariateMomentSet> mPtr
    (
        new univariateMomentSet(nMoments, support, SMALL, 0.0)
    );

    mPtr().update(weights, abscissae);

    return mPtr;
}


//- Take a moment set through the auxiliary quantities and back, and require
//  that it comes back unchanged
void checkRoundTrip
(
    const scalarList& weights,
    const scalarList& abscissae,
    const label nMoments,
    const supportType& support,
    const scalar tolerance,
    const string& what
)
{
    autoPtr<univariateMomentSet> mPtr
    (
        momentSetFromQuadrature(weights, abscissae, nMoments, support)
    );

    univariateMomentSet& m = mPtr();

    const label nRealizable = m.nRealizableMoments(false);

    check
    (
        nRealizable == nMoments,
        what + ": the moment set built from the quadrature is realizable"
    );

    scalarList expected(nMoments);

    forAll(expected, mi)
    {
        expected[mi] = m[mi];
    }

    const scalar m0 = expected[0];

    scalarList recovered(nMoments, Zero);

    if (support == supportType::RPlus)
    {
        univariateAdvection::zeta::zetaToMoments
        (
            m.zetas(), recovered, nMoments, m0
        );
    }
    else
    {
        univariateAdvection::zeta::canonicalMomentsToMoments
        (
            m.canonicalMoments(), recovered, nMoments, m0
        );
    }

    checkMoments(recovered, expected, tolerance, what);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- A measure over R+ recovered from its zeta chain
void testZetaRoundTripOnRPlus()
{
    Info<< "\nMoments over R+ recovered from the zeta chain" << endl;

    // Two nodes, a unit zero-order moment
    checkRoundTrip
    (
        {0.5, 0.5}, {1.0, 2.0}, 4, supportType::RPlus, 1.0e-12,
        "two equally weighted nodes, four moments"
    );

    // Two nodes, unequal weights, a zero-order moment away from one
    checkRoundTrip
    (
        {0.3, 1.7}, {0.5, 3.0}, 4, supportType::RPlus, 1.0e-12,
        "two unequally weighted nodes, four moments"
    );

    // Three nodes, six moments. This is the size at which prod[] is indexed
    // one past the end of the zeta chain in zetaToMoments
    checkRoundTrip
    (
        {0.2, 0.5, 0.3}, {0.5, 1.5, 4.0}, 6, supportType::RPlus, 1.0e-11,
        "three nodes, six moments"
    );

    // Four nodes, eight moments
    checkRoundTrip
    (
        {0.1, 0.4, 0.3, 0.2}, {0.25, 1.0, 2.5, 6.0}, 8,
        supportType::RPlus, 1.0e-9,
        "four nodes, eight moments"
    );

    // Abscissae spanning several decades, which stresses the zeta chain
    checkRoundTrip
    (
        {0.25, 0.5, 0.25}, {1.0e-3, 1.0, 1.0e2}, 6,
        supportType::RPlus, 1.0e-9,
        "three nodes spanning five decades, six moments"
    );
}


//- A measure over [0, 1] recovered from its canonical moments
void testCanonicalMomentRoundTripOnUnitInterval()
{
    Info<< "\nMoments over [0, 1] recovered from the canonical moments" << endl;

    checkRoundTrip
    (
        {0.5, 0.5}, {0.25, 0.75}, 4, supportType::ZeroOne, 1.0e-12,
        "two equally weighted nodes, four moments"
    );

    checkRoundTrip
    (
        {0.7, 0.3}, {0.1, 0.9}, 4, supportType::ZeroOne, 1.0e-12,
        "two unequally weighted nodes, four moments"
    );

    checkRoundTrip
    (
        {0.2, 0.5, 0.3}, {0.05, 0.4, 0.95}, 6,
        supportType::ZeroOne, 1.0e-11,
        "three nodes, six moments"
    );

    // Nodes crowded against the ends of the interval, where the canonical
    // moments approach their bounds
    checkRoundTrip
    (
        {0.5, 0.5}, {1.0e-4, 1.0 - 1.0e-4}, 4,
        supportType::ZeroOne, 1.0e-10,
        "two nodes at the ends of the interval, four moments"
    );
}


//- The zero-order moment is carried separately by the scheme, so the map has
//  to be linear in it
void testZeroOrderMomentScaling()
{
    Info<< "\nScaling with the zero-order moment" << endl;

    const scalarList weights({0.2, 0.5, 0.3});
    const scalarList abscissae({0.5, 1.5, 4.0});
    const label nMoments = 6;

    autoPtr<univariateMomentSet> mPtr
    (
        momentSetFromQuadrature
        (
            weights, abscissae, nMoments, supportType::RPlus
        )
    );

    univariateMomentSet& m = mPtr();
    m.nRealizableMoments(false);

    scalarList unit(nMoments, Zero);
    univariateAdvection::zeta::zetaToMoments(m.zetas(), unit, nMoments, 1.0);

    const scalar m0 = 17.0;

    scalarList scaled(nMoments, Zero);
    univariateAdvection::zeta::zetaToMoments(m.zetas(), scaled, nMoments, m0);

    scalarList expected(nMoments);

    forAll(expected, mi)
    {
        expected[mi] = m0*unit[mi];
    }

    checkMoments
    (
        scaled, expected, 1.0e-12,
        "the reconstructed moments are linear in the zero-order moment"
    );
}


//- A distribution carried by a single node is the degenerate case the scheme
//  meets in a cell holding one size, and the map has to reproduce it exactly
void testSingleNode()
{
    Info<< "\nA distribution carried by a single node" << endl;

    const scalar x = 2.5;
    const label nMoments = 6;

    scalarList expected(nMoments);

    forAll(expected, mi)
    {
        expected[mi] = pow(x, mi);
    }

    autoPtr<univariateMomentSet> mPtr
    (
        momentSetFromQuadrature
        (
            {1.0}, {x}, nMoments, supportType::RPlus
        )
    );

    univariateMomentSet& m = mPtr();
    m.nRealizableMoments(false);

    scalarList recovered(nMoments, Zero);
    univariateAdvection::zeta::zetaToMoments
    (
        m.zetas(), recovered, nMoments, 1.0
    );

    checkMoments
    (
        recovered, expected, 1.0e-12,
        "a single node reproduces the moments of a Dirac measure"
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    testZetaRoundTripOnRPlus();
    testCanonicalMomentRoundTripOnUnitInterval();
    testZeroOrderMomentScaling();
    testSingleNode();

    Info<< "\n" << nTested << " checks passed.\n" << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
