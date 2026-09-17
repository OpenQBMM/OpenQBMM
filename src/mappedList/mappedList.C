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

\*---------------------------------------------------------------------------*/

#include "mappedList.H"


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template <class mappedType>
Foam::word
Foam::mappedList<mappedType>::listToWord(const labelList& list)
{
    return mappedListKey::ordersToWord(list);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class mappedType>
void Foam::mappedList<mappedType>::setIndexes(const labelListList& indexes)
{
    if (indexes.size() != this->size())
    {
        FatalErrorInFunction
            << "Size mismatch: "
            << endl
            << "  List size = " << this->size()
            << endl
            << "  Indexes list size = " << indexes.size()
            << exit(FatalError);
    }

    index_ = mappedListKey(indexes);
}


template <class mappedType>
void Foam::mappedList<mappedType>::notFound(const labelList& list) const
{
    FatalErrorInFunction
        << "No element of orders " << list << " in a mapped list of "
        << index_.nDimensions() << " dimensions." << nl
        << exit(FatalError);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const label size,
    const labelListList& indexes
)
:
    List<mappedType>(size),
    index_()
{
    setIndexes(indexes);
}

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const label size,
    const labelListList& indexes,
    const mappedType& initValue
)
:
    List<mappedType>(size, initValue),
    index_()
{
    setIndexes(indexes);
}

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const List<mappedType>& initList,
    const labelListList& indexes
)
:
    List<mappedType>(initList),
    index_()
{
    setIndexes(indexes);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class mappedType>
void Foam::mappedList<mappedType>::setSize
(
    const label newSize,
    const labelListList& newIndexes
)
{
    Foam::List<mappedType>::setSize(newSize);
    setIndexes(newIndexes);
}

template <class mappedType>
void Foam::mappedList<mappedType>::resize
(
    const label newSize,
    const labelListList& newIndexes
)
{
    (*this).setSize(newSize, newIndexes);
}

template <class mappedType>
bool Foam::mappedList<mappedType>::found(const labelList& list) const
{
    // Access takes zeros past the dimensions of the list as the element
    // without them; found is asked whether the orders, as given, belong
    // to the list
    return list.size() <= index_.nDimensions() && index_.find(list) >= 0;
}

template <class mappedType>
template <typename ...ArgsT>
bool Foam::mappedList<mappedType>::found(ArgsT...args) const
{
    return
        label(sizeof...(args)) <= index_.nDimensions()
     && index_.find({args...}) >= 0;
}

// ************************************************************************* //
