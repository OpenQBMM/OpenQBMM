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
    Foam::mappedList

Description
    A list with additional mapping for named access to elements.

    Every element is identified by a list of orders, or indexes, given when
    the list is built, and can be reached by them as well as by its position:

        Orders  |  Position
         0 0 0         0
         1 0 0         1
         0 1 0         2
         0 0 1         3
         2 0 0         4
         2 1 0         5

    An element is found by its orders exactly, whatever their values and
    however many dimensions they have (see mappedListKey). Orders followed by
    zeros are the same element as the orders without them, so moments(1)
    reaches the element of orders (1 0 0) of a list of three dimensions.

    Provides O(1) access by position and by orders.

SourceFiles
    mappedListI.H
    mappedList.C

\*---------------------------------------------------------------------------*/

#ifndef mappedList_H
#define mappedList_H

#include "scalar.H"
#include "vector.H"
#include "labelList.H"
#include "mappedListKey.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mappedList Declaration
\*---------------------------------------------------------------------------*/

template <class mappedType>
class mappedList : public List<mappedType>
{
    // Private data

        //- The orders of every element and the table that finds them
        mappedListKey index_;


    // Private Member Functions

        //- Set the orders of every element, checking the size of the list
        void setIndexes(const labelListList& indexes);

        //- Stop on orders that belong to no element
        void notFound(const labelList& list) const;


public:

    // Constructors

        //- Construct from size and the orders of each element
        mappedList
        (
            const label size,
            const labelListList& indexes
        );

        //- Construct from size and the orders of each element, with an
        //  initial value for all elements
        mappedList
        (
            const label size,
            const labelListList& indexes,
            const mappedType& initValue
        );

        //- Construct from List<mappedType> and the orders of each element
        mappedList
        (
            const List<mappedType>& initList,
            const labelListList& indexes
        );


    // Member Functions

        // Static member functions

            //- Convert a list of labels to a word. {1, 2, 3} -> 1_2_3, {3} -> 3
            static word listToWord(const labelList& list);


        // Modifiers

            //- Set size and rebuild map with new indexes
            void setSize
            (
                const label newSize,
                const labelListList& newIndexes
            );

            //- Resize and rebuild map with new indexes
            void resize
            (
                const label newSize,
                const labelListList& newIndexes
            );


        // Access

            //- Is label list index used
            bool found(const labelList& list) const;

            //- Is label list index used
            template <typename ...ArgsT>
            bool found(ArgsT ...args) const;

            //- Non-const access to a moment given the list of orders
            template <typename ...ArgsT>
            inline mappedType& operator()(ArgsT ...args);

            //- Const access to a moment given the list of orders
            template <typename ...ArgsT>
            inline const mappedType& operator()(ArgsT ...args) const;

            //- Constant access given a label list
            inline const mappedType& operator()(const labelList& list) const;

            //- Non-constant access given a label list
            inline mappedType& operator()(const labelList& list);

            //- The orders of each element, by position, padded with zeros
            //  to the number of dimensions
            inline labelListList indexes() const;

            //- Number of dimensions of the orders
            inline label nDimensions() const noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mappedListI.H"
#   include "mappedList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
