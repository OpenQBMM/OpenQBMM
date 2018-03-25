/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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
    Foam::momentSet

Description
    Stores a generic moment set with the corresponding list of orders.

SourceFiles
    momentSet.C
    momentSetI.H

\*---------------------------------------------------------------------------*/

#ifndef momentSet_H
#define momentSet_H

#include "labelList.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class momentSet Declaration
\*---------------------------------------------------------------------------*/

class momentSet
:
    public scalarList
{
protected:
    // Protected data

        //- Number of moments
        label nMoments_;

        //- List of moment orders
        labelListList momentOrders_;

        //- Type of support of the measure associated to the moment set
        word support_;


public:

    // Constructors

        //- Construct from number of moments and initial value
        momentSet
        (
            const label nMoments,
            const labelListList& momentOrders,
            const word& support,
            const scalar initValue = 0
        );

        //- Construct from scalarList
        momentSet
        (
            const scalarList& m,
            const labelListList& momentOrders,
            const word& support
        );

    //- Destructor
    virtual ~momentSet();


    // Member Functions

        //- Return moment orders
        inline const labelListList& momentOrders() const;

        //- Return the number of moments
        inline label nMoments() const;

        //- Return support
        inline const word support() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "momentSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
