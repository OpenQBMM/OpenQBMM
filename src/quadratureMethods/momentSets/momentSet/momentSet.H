/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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
#include "mappedLists.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class momentSet Declaration
\*---------------------------------------------------------------------------*/

class momentSet : public mappedList<scalar>
{
protected:

    // Protected data

        //- Number of moments
        label nMoments_;

        //- Number of dimensions of the the NDF associated to the moment set
        label nDimensions_;

        //- List of moment orders
        labelListList momentOrders_;

        //- List of supports of the measure associated to the moment set
        wordList supports_;

        //- Minimum value for the zero-order moment.
        // Below this value, the moment is considered negligible.
        scalar smallM0_;

        //- Minimum value for zeta_k.
        // Below this value, zeta_k is considered negligible.
        scalar smallZeta_;


    // Protected Member Functions

        //- Validate list of supports
        void validateSupports() const;

        //- Validate compatibility between number of dimensions and moment orders
        void validateMomentOrders() const;

        //- Validate
        void validate() const;

public:

    // Constructors

        //- Construct from number of moments and initial value
        momentSet
        (
            const label nMoments,
            const label nDimensions,
            const labelListList& momentOrders,
            const wordList& supports,
            const scalar smallM0,
            const scalar smallZeta,
            const scalar initialValue = 0
        );

        //- Construct from scalarList
        momentSet
        (
            const scalarList& m,
            const label nDimensions,
            const labelListList& momentOrders,
            const wordList& supports,
            const scalar smallM0,
            const scalar smallZeta
        );

    //- Destructor
    virtual ~momentSet();


    // Member Functions

        //- Return moment orders
        inline const labelListList& momentOrders() const;

        //- Return the number of dimensions of the NDF
        inline label nDimensions() const;

        //- Return the number of moments
        inline label nMoments() const;

        //- Resizes moment set
        void resize(const label newSize);

        //- Set size of moment set
        void setSize(const label newSize);

        //- Return support
        inline const wordList& supports() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "momentSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
