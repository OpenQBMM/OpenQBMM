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
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    word listWord;

    forAll(list, dimi)
    {
        listWord += Foam::name(list[dimi]);
    }

    return listWord;
}

template <class mappedType>
Foam::label
Foam::mappedList<mappedType>::listToLabel
(
    const labelList& list,
    const label nDimensions
)
{
    label listLabel = 0;
    label size = max(nDimensions, list.size());

    forAll(list, dimi)
    {
        listLabel += list[dimi]*pow(scalar(10), size - dimi - 1);
    }

    return listLabel;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const label size,
    const labelListList& indexes
)
:
    List<mappedType>(size),
    map_(size),
    nDimensions_(0)
{
    forAll(indexes, i)
    {
        nDimensions_ = max(nDimensions_, indexes[i].size());
    }

    forAll(*this, elemi)
    {
        map_.insert
        (
            listToLabel(indexes[elemi], nDimensions_),
            elemi
        );
    }
}

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const label size,
    const labelListList& indexes,
    const mappedType& initValue
)
:
    List<mappedType>(size, initValue),
    map_(size),
    nDimensions_(0)
{
    forAll(indexes, i)
    {
        nDimensions_ = max(nDimensions_, indexes[i].size());
    }

    forAll(*this, elemi)
    {
        map_.insert
        (
            listToLabel(indexes[elemi], nDimensions_),
            elemi
        );
    }
}

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const label size,
    const Map<label>& map,
    const mappedType& initValue
)
:
    List<mappedType>(size, initValue),
    map_(map),
    nDimensions_(0)
{
    forAllConstIter(Map<label>, map_, iter)
    {
        label key = iter.key();
        label nD = 0;

        while (key)
        {
            key /= 10;
            nD++;
        }

        nDimensions_ = max(nDimensions_, nD);
    }
}

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const List<mappedType>& initList,
    const labelListList& indexes
)
:
    List<mappedType>(initList),
    map_(initList.size()),
    nDimensions_(0)
{
    forAll(indexes, i)
    {
        nDimensions_ = max(nDimensions_, indexes[i].size());
    }

    forAll(*this, elemi)
    {
        map_.insert
        (
            listToLabel(indexes[elemi], nDimensions_),
            elemi
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class mappedType>
Foam::mappedList<mappedType>::~mappedList()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class mappedType>
Foam::label Foam::mappedList<mappedType>::calcMapIndex
(
    std::initializer_list<Foam::label> indexes
) const
{
    label mapIndex = 0;

    if (indexes.size() > 0)
    {
        for
        (
            std::initializer_list<label>::iterator iter = indexes.begin();
            iter < indexes.end();
            iter++
        )
        {
            label argIndex = std::distance(indexes.begin(), iter);
            mapIndex += (*iter)*pow(scalar(10), nDimensions_ - argIndex - 1);
        }
    }

    return mapIndex;
}

template <class mappedType>
void Foam::mappedList<mappedType>::setSize(const label newSize)
{
    Foam::List<mappedType>::setSize(newSize);
    map_.resize(newSize);
}

template <class mappedType>
void Foam::mappedList<mappedType>::resize(const label newSize)
{
    (*this).setSize(newSize);
}

template <class mappedType>
bool Foam::mappedList<mappedType>::found(const labelList& list) const
{
    if (list.size() > nDimensions_)
    {
        return false;
    }

    forAllConstIter(Map<label>, map_, iter)
    {
        label key = iter.key();

        if (key == listToLabel(list, nDimensions_))
        {
            return true;
        }
    }

    return false;
}

template <class mappedType>
template <typename ...ArgsT>
bool Foam::mappedList<mappedType>::found(ArgsT...args) const
{
    if (label(std::initializer_list<Foam::label>({args...}).size()) > nDimensions_)
    {
        return false;
    }

    forAllConstIter(Map<label>, map_, iter)
    {
        label key = iter.key();

        if (key == calcMapIndex({args...}))
        {
            return true;
        }
    }

    return false;
}

// ************************************************************************* //
