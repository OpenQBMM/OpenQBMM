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

Class
    Foam::realizableOdeSolver

Description
    Solve a univariate PDF transport equation with the extended quadrature
    method of moments.

SourceFiles
    realizableOdeSolver.C

\*---------------------------------------------------------------------------*/

#ifndef realizableOdeSolver_H
#define realizableOdeSolver_H

#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class realizableOdeSolver Declaration
\*---------------------------------------------------------------------------*/

template<class momentType, class nodeType>
class realizableOdeSolver
{

    typedef momentFieldSet<momentType, nodeType> momentFieldSetType;
    typedef quadratureApproximation<momentFieldSetType, nodeType> quadratureType;

    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Ode coefficient dictionary
        const dictionary& odeDict_;

        //- Absolute tolerance for realizable ODE solver
        scalar ATol_;

        //- Relative tolerance for realizable ODE solver
        scalar RTol_;

        //- Time-step change factor
        scalar fac_;

        //- Minimum time-step change factor
        scalar facMin_;

        //- Maximum time-step change factor
        scalar facMax_;

        //- Minimum local step in realizable ODE solver
        scalar minLocalDt_;

        //- Stored time step values
        scalarField localDt_;


protected:

    // Protected data

        //- Are the moments updated due to sources
        Switch solveSources_;

        //- Is the realizable ode solver used
        Switch solveOde_;

    // Protected member functions

        //- Calculate source term for moment equation in the specified cell
        virtual void updateCellMomentSource(const label celli) = 0;

        //- Calculate source term for moment equation in the specified cell
        virtual scalar cellMomentSource
        (
            const label momentOrder,
            const label celli,
            const mappedPtrList<nodeType>& nodes,
            const label enviroment
        ) = 0;


public:

    // Constructors

        //- Construct from components
        realizableOdeSolver
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~realizableOdeSolver();


    // Member Functions

        //- Solve univariate PDF transport equation
        void solve(quadratureType& quadrature, const label enviroment);

        //- Read coefficients from dictionary
        void read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if NoRepository
    #include "realizableOdeSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif