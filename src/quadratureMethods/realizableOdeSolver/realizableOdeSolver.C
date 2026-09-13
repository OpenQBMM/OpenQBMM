/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/

#include "realizableOdeSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class momentType, class nodeType>
Foam::realizableOdeSolver<momentType, nodeType>::realizableOdeSolver
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    ATol_(readScalar(dict.subDict("odeCoeffs").lookup("ATol"))),
    RTol_(readScalar(dict.subDict("odeCoeffs").lookup("RTol"))),
    fac_(readScalar(dict.subDict("odeCoeffs").lookup("fac"))),
    facMin_(readScalar(dict.subDict("odeCoeffs").lookup("facMin"))),
    facMax_(readScalar(dict.subDict("odeCoeffs").lookup("facMax"))),
    minLocalDt_(readScalar(dict.subDict("odeCoeffs").lookup("minLocalDt"))),
    localDt_
    (
        IOobject
        (
            "realizableOde:localDt",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        mesh.time().deltaT()
    ),
    localDtAdjustments_(0),
    solveSources_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("solveSources", true)
    ),
    solveOde_
    (
        dict.subDict("odeCoeffs").lookupOrDefault("solveOde", true)
    )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class momentType, class nodeType>
Foam::realizableOdeSolver<momentType, nodeType>::~realizableOdeSolver()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class momentType, class nodeType>
void Foam::realizableOdeSolver<momentType, nodeType>::solve
(
    quadratureType& quadrature,
    const label enviroment
)
{
    if (!solveSources_)
    {
        return;
    }

    momentFieldSetType& moments(quadrature.moments());
    label nMoments = quadrature.nMoments();
    scalar globalDt = mesh_.time().deltaT().value();
    const labelListList& momentOrders = quadrature.momentOrders();

    //- Use Euler explicit to update moments due to sources
    if (!solveOde_)
    {
        forAll(moments[0], celli)
        {
            updateCellMomentSource(celli);

            forAll(moments, mi)
            {
                const labelList& order = momentOrders[mi];

                moments[mi][celli] +=
                    globalDt
                   *cellMomentSource
                    (
                        order,
                        celli,
                        quadrature,
                        enviroment
                    );
            }

            quadrature.updateLocalQuadrature(celli, true);
            quadrature.updateLocalMoments(celli);
        }

        forAll(moments, mi)
        {
            moments[mi].correctBoundaryConditions();
        }

        quadrature.updateBoundaryQuadrature();

        return;
    }

    Info << "Solving source terms in realizable ODE solver." << endl;

    forAll(moments[0], celli)
    {
        //Info << "OLD MOMENTS" << moments << endl;

        // Storing old moments to recover from failed step
        quadrature.updateLocalQuadrature(celli);
        quadrature.updateLocalMoments(celli);

        scalarList oldMoments(nMoments, Zero);

        forAll(oldMoments, mi)
        {
            oldMoments[mi] = moments[mi][celli];
        }

        //Info << "Old moments: " << oldMoments << endl;

        //- Local time
        scalar localT(0);

        scalar localDt = min(localDt_[celli], globalDt);

        // Initialize RK parameters
        scalarList k1(nMoments, Zero);
        scalarList k2(nMoments, Zero);
        scalarList k3(nMoments, Zero);

        // Flag to indicate if the time step is complete
        bool timeComplete = false;

        // Check realizability of intermediate moment sets
        bool realizableUpdate1 = false;
        bool realizableUpdate2 = false;
        bool realizableUpdate3 = false;

        scalarList diff23(nMoments, Zero);

        while (!timeComplete)
        {
            do
            {
                // First intermediate update
                bool nullSource =  true;

                updateCellMomentSource(celli);

                forAll(k1, mi)
                {
                    const labelList& order = momentOrders[mi];

                    k1[mi] =
                        localDt*cellMomentSource
                        (
                            order,
                            celli,
                            quadrature,
                            enviroment
                        );

                    moments[mi][celli] = oldMoments[mi] + k1[mi];

                    // A source is null when it moves no moment by more than
                    // a vanishing part of the tolerance of that moment. It
                    // used to be null below SMALL, an absolute number: a
                    // cell whose local step had shrunk to round-off then had
                    // every source declared null, and stayed frozen.
                    if
                    (
                        mag(k1[mi])
                      > SMALL*(ATol_ + RTol_*mag(oldMoments[mi]))
                    )
                    {
                        nullSource = false;
                    }
                }

                // Nothing acts on the cell, so it stays where it is for the
                // whole of the step. The moments are put back rather than
                // left with the negligible k1, and the quadrature with them.
                if (nullSource)
                {
                    forAll(oldMoments, mi)
                    {
                        moments[mi][celli] = oldMoments[mi];
                    }

                    quadrature.updateLocalQuadrature(celli);

                    localT = globalDt;
                    timeComplete = true;
                    break;
                }

                realizableUpdate1 =
                        quadrature.updateLocalQuadrature(celli, false);

                quadrature.updateLocalMoments(celli);

                // Second moment update
                updateCellMomentSource(celli);

                forAll(k2, mi)
                {
                    const labelList& order = momentOrders[mi];

                    k2[mi] =
                        localDt*cellMomentSource
                        (
                            order,
                            celli,
                            quadrature,
                            enviroment
                        );

                    moments[mi][celli] = oldMoments[mi] + (k1[mi] + k2[mi])/4.0;
                }

                realizableUpdate2 =
                    quadrature.updateLocalQuadrature(celli, false);

                quadrature.updateLocalMoments(celli);

                // Third moment update
                updateCellMomentSource(celli);

                forAll(k3, mi)
                {
                    const labelList& order = momentOrders[mi];

                    k3[mi] =
                        localDt*cellMomentSource
                        (
                            order,
                            celli,
                            quadrature,
                            enviroment
                        );

                    moments[mi][celli] =
                        oldMoments[mi] + (k1[mi] + k2[mi] + 4.0*k3[mi])/6.0;

                    diff23[mi] = (2.0*k3[mi] - k1[mi] - k2[mi])/3.0;
                }

                realizableUpdate3 =
                    quadrature.updateLocalQuadrature(celli, false);

                quadrature.updateLocalMoments(celli);

                if
                (
                    !realizableUpdate1
                 || !realizableUpdate2
                 || !realizableUpdate3
                )
                {
                    // Avoid spamming the terminal when not realizable
                    if (localDtAdjustments_ == 0)
                    {
                        Info << "Not realizable, adjusting local timestep."
                             << nl
                             << "This may take a while." << endl;
                    }

                    localDtAdjustments_++;

                    forAll(oldMoments, mi)
                    {
                        moments[mi][celli] = oldMoments[mi];
                    }

                    // Updating local quadrature with old moments
                    quadrature.updateLocalQuadrature(celli);

                    localDt /= 2.0;

                    if (localDt < minLocalDt_)
                    {
                        FatalErrorInFunction
                            << "Reached minimum local step in realizable ODE"
                            << nl
                            << "    solver. Cannot ensure realizability."
                            << nl
                            << abort(FatalError);
                    }
                }
            }
            while
            (
                !realizableUpdate1
             || !realizableUpdate2
             || !realizableUpdate3
            );

            if (timeComplete)
            {
                break;
            }

            // Initialize error and change
            scalar error(0);
            scalar maxChange(0);

            for (label mi = 0; mi < nMoments; mi++)
            {
                // Calculate the scaling factor
                scalar scalei =
                    ATol_
                  + max
                    (
                        mag(moments[mi][celli]), mag(oldMoments[mi])
                    )*RTol_;

                // Update the error
                error += sqr(diff23[mi]/scalei);

                // Update the largest change in the moments, measured
                // against the same scale as the error. Taken as it comes,
                // the change is a moment difference, so what counts as
                // small depends on how large the moments of the case are:
                // a cell whose moments are of order 1e-12 was leaving the
                // integration on a change that is a part in a thousand of
                // them, while one of order 1e6 could never reach the
                // threshold at all. The scale is what ATol and RTol are
                // for, and the error beside it already uses it.
                maxChange =
                    max
                    (
                        maxChange,
                        mag(moments[mi][celli] - oldMoments[mi])/scalei
                    );
            }

            error = sqrt(error/nMoments);

            // A substep over which no moment changed, to a vanishing part of
            // its tolerance, leaves a cell that is stationary: it is done
            // for the whole of the step. This used to be taken as complete
            // as soon as either the change or the error estimate was
            // vanishing, without the time integrated being advanced, so the
            // rest of the step was dropped; and a substep the controller
            // had already shrunk to round-off met it at once, every step.
            if (maxChange < SMALL)
            {
                localT = globalDt;
                timeComplete = true;

                if (error >= 1)
                {
                    WarningInFunction
                        << "The maximum change in moments is small, "
                        << "but error is not.\n"
                        << nl
                        << "Error: " << error << nl
                        << "Max. scaled change: " << maxChange << nl
                        << nl
                        << "\nThis may indicate a problem with the "
                        << "realizable ODE solver." << endl;
                }

                break;
            }
            else if (error < 1)
            {
                localT += localDt;

                localDt *= min(facMax_, max(facMin_, fac_/pow(error, 1.0/3.0)));

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }

                // What is left of the global step. A remainder that is a
                // vanishing part of it is round-off accumulated over the
                // substeps, and is taken as the end of the step.
                const scalar remaining = globalDt - localT;

                if (remaining <= ROOTSMALL*globalDt)
                {
                    localT = globalDt;
                    timeComplete = true;
                    break;
                }

                localDt = min(remaining, localDt);
                localDt_[celli] = localDt;
            }
            else
            {
                localDt *=
                    min(scalar(1), max(facMin_, fac_/pow(error, 1.0/3.0)));

                forAll(oldMoments, mi)
                {
                    moments[mi][celli] = oldMoments[mi];
                }

                // Updating local quadrature with old moments
                quadrature.updateLocalQuadrature(celli);
            }
        }
    }

    forAll(moments, mi)
    {
        moments[mi].correctBoundaryConditions();
    }

    quadrature.updateBoundaryQuadrature();
}


template<class momentType, class nodeType>
void Foam::realizableOdeSolver<momentType, nodeType>
::read(const dictionary& dict)
{
    const dictionary& odeDict = dict.subDict("odeCoeffs");
    solveSources_ = odeDict.lookupOrDefault<Switch>("solveSources", true);
    solveOde_ = odeDict.lookupOrDefault<Switch>("solveOde", true);

    (odeDict.lookup("ATol")) >> ATol_;
    (odeDict.lookup("RTol")) >> RTol_;
    (odeDict.lookup("fac")) >> fac_;
    (odeDict.lookup("facMin")) >> facMin_;
    (odeDict.lookup("facMax")) >> facMax_;
    (odeDict.lookup("minLocalDt")) >> minLocalDt_;
}


// ************************************************************************* //
