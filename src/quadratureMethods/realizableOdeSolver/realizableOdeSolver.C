/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
     \\/     M anipulation  |
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
    odeDict_(dict.subDict("odeCoeffs")),
    ATol_(readScalar(odeDict_.lookup("ATol"))),
    RTol_(readScalar(odeDict_.lookup("RTol"))),
    fac_(readScalar(odeDict_.lookup("fac"))),
    facMin_(readScalar(odeDict_.lookup("facMin"))),
    facMax_(readScalar(odeDict_.lookup("facMax"))),
    minLocalDt_(readScalar(odeDict_.lookup("minLocalDt"))),
    localDt_(mesh.nCells(), mesh.time().deltaTValue()/10.0),
    solveSources_
    (
        odeDict_.lookupOrDefault("solveSources", true)
    ),
    solveOde_
    (
        odeDict_.lookupOrDefault("solveOde", true)
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
    const mappedPtrList<nodeType>& nodes(quadrature.nodes());
    label nMoments = quadrature.nMoments();
    scalar globalDt = mesh_.time().deltaT().value();

    if (!solveOde_)
    {
        forAll(moments[0], celli)
        {
            updateCellMomentSource(celli);
            forAll(moments, mi)
            {
                moments[mi][celli] +=
                    globalDt*cellMomentSource(mi, celli, nodes, enviroment);
            }

            quadrature.updateLocalQuadrature(celli, true);
            quadrature.updateLocalMoments(celli);
        }
        return;
    }

    Info << "Solving source terms in realizable ODE solver." << endl;

    forAll(moments[0], celli)
    {
        // Storing old moments to recover from failed step

        scalarList oldMoments(nMoments, 0.0);

        forAll(oldMoments, mi)
        {
            oldMoments[mi] = moments[mi][celli];
        }

        //- Local time
        scalar localT = 0.0;

        // Initialize the local step
        scalar localDt = localDt_[celli];

        // Initialize RK parameters
        scalarList k1(nMoments, 0.0);
        scalarList k2(nMoments, 0.0);
        scalarList k3(nMoments, 0.0);

        // Flag to indicate if the time step is complete
        bool timeComplete = false;

        // Check realizability of intermediate moment sets
        bool realizableUpdate1 = false;
        bool realizableUpdate2 = false;
        bool realizableUpdate3 = false;

        scalarList momentsSecondStep(nMoments, 0.0);

        while (!timeComplete)
        {
            do
            {
                // First intermediate update
                updateCellMomentSource(celli);
                forAll(oldMoments, mi)
                {
                    k1[mi] =
                        localDt*cellMomentSource(mi, celli, nodes, enviroment);
                    moments[mi][celli] = oldMoments[mi] + k1[mi];
                }

                realizableUpdate1 =
                        quadrature.updateLocalQuadrature(celli, false);

                quadrature.updateLocalMoments(celli);

                // Second moment update
                updateCellMomentSource(celli);
                forAll(oldMoments, mi)
                {
                    k2[mi] =
                        localDt*cellMomentSource(mi, celli, nodes, enviroment);
                    moments[mi][celli] = oldMoments[mi] + (k1[mi] + k2[mi])/4.0;

                    momentsSecondStep[mi] = moments[mi][celli];
                }

                realizableUpdate2 =
                        quadrature.updateLocalQuadrature(celli, false);

                quadrature.updateLocalMoments(celli);

                // Third moment update
                updateCellMomentSource(celli);
                forAll(oldMoments, mi)
                {
                    k3[mi] =
                        localDt*cellMomentSource(mi, celli, nodes, enviroment);
                    moments[mi][celli] =
                        oldMoments[mi] + (k1[mi] + k2[mi] + 4.0*k3[mi])/6.0;
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
                    Info << "Not realizable" << endl;

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
                            << "    solver. Cannot ensure realizability." << nl
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

            scalar error = 0.0;

            for (label mi = 0; mi < nMoments; mi++)
            {
                scalar scalei =
                        ATol_
                    + max
                        (
                            mag(momentsSecondStep[mi]), mag(oldMoments[mi])
                        )*RTol_;

                error +=
                        sqr
                        (
                            (momentsSecondStep[mi] - moments[mi][celli])/scalei
                        );
            }

            error = max(sqrt(error/nMoments), small);

            if (error == 0)
            {
                timeComplete = true;
                localT = 0.0;
                break;
            }
            else if (error < 1)
            {
                localDt *= min(facMax_, max(facMin_, fac_/pow(error, 1.0/3.0)));

                scalar maxLocalDt = max(globalDt - localT, 0.0);
                localDt = min(maxLocalDt, localDt);

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }

                if (localDt == 0.0)
                {
                    timeComplete = true;
                    localT = 0.0;
                    break;
                }

                localDt_[celli] = localDt;
                localT += localDt;
            }
            else
            {
                localDt *= min(1.0, max(facMin_, fac_/pow(error, 1.0/3.0)));

                forAll(oldMoments, mi)
                {
                    moments[mi][celli] = oldMoments[mi];
                }

                // Updating local quadrature with old moments
                quadrature.updateLocalQuadrature(celli);
            }
        }
    }
}


template<class momentType, class nodeType>
void Foam::realizableOdeSolver<momentType, nodeType>::read()
{
    odeDict_.lookupOrDefault("solveSources", true) >> solveSources_;
    odeDict_.lookupOrDefault("solveOde", true) >> solveOde_;

    (odeDict_.lookup("ATol")) >> ATol_;
    (odeDict_.lookup("RTol")) >> RTol_;
    (odeDict_.lookup("fac")) >> fac_;
    (odeDict_.lookup("facMin")) >> facMin_;
    (odeDict_.lookup("facMax")) >> facMax_;
    (odeDict_.lookup("minLocalDt")) >> minLocalDt_;
}


// ************************************************************************* //
