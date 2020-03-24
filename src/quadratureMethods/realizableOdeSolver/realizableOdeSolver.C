/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2020 Alberto Passalacqua
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
        // Storing old moments to recover from failed step
        quadrature.updateLocalQuadrature(celli);
        quadrature.updateLocalMoments(celli);

        scalarList oldMoments(nMoments, Zero);
        forAll(oldMoments, mi)
        {
            oldMoments[mi] = moments[mi][celli];
        }

        //- Local time
        scalar localT(0);

        // Initialize the local step
        scalar localDt = localDt_[celli];

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
        label nItt = 0;

        while (!timeComplete)
        {
            do
            {
                nItt++;

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

                    if (mag(k1[mi]) > SMALL)
                    {
                        nullSource = false;
                    }
                }

                realizableUpdate1 =
                        quadrature.updateLocalQuadrature(celli, false);

               quadrature.updateLocalMoments(celli);

               if (nullSource)
               {
                   break;
               }

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

            scalar error(0);

            for (label mi = 0; mi < nMoments; mi++)
            {
                scalar scalei =
                    ATol_
                  + max
                    (
                        mag(moments[mi][celli]), mag(oldMoments[mi])
                    )*RTol_;

                error += sqr(diff23[mi]/scalei);
            }

            error = sqrt(error/nMoments);

            if (error < SMALL)
            {
                timeComplete = true;
                localT = Zero;
                break;
            }
            else if (error < 1)
            {
                localT += localDt;
                localDt *= min(facMax_, max(facMin_, fac_/pow(error, 1.0/3.0)));

                scalar maxLocalDt = max(globalDt - localT, scalar(0));
                localDt = min(maxLocalDt, localDt);

                forAll(oldMoments, mi)
                {
                    oldMoments[mi] = moments[mi][celli];
                }

                if (localDt == 0.0)
                {
                    timeComplete = true;
                    localT = Zero;
                    break;
                }

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
