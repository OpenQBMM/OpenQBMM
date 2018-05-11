/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
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

#include "mappedList.H"


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template <class mappedType>
Foam::word
Foam::mappedList<mappedType>::listToWord(const labelList& lst)
{
    word w;

    forAll(lst, dimi)
    {
        w += Foam::name(lst[dimi]);
    }

    return w;
}

template <class mappedType>
Foam::label
Foam::mappedList<mappedType>::listToLabel(const labelList& lst)
{
    label l = 0;

    forAll(lst, dimi)
    {
        l += lst[dimi]*pow(10, lst.size() - dimi - 1);
    }

    return l;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const label size,
    const labelListList& indexes
)
:
    List<mappedType>(size),
    map_(size)
{
    forAll(*this, elemi)
    {
        map_.insert
        (
            listToLabel(indexes[elemi]),
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
    map_(size)
{
    forAll(*this, elemi)
    {
        map_.insert
        (
            listToLabel(indexes[elemi]),
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
    map_(map)
{}

template <class mappedType> Foam::mappedList<mappedType>::mappedList
(
    const List<mappedType>& initList,
    const labelListList& indexes
)
:
    List<mappedType>(initList),
    map_(initList.size())
{
    forAll(*this, elemi)
    {
        map_.insert
        (
            listToLabel(indexes[elemi]),
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
    label argSize = indexes.size();
    label mapIndex = 0;

    if (argSize > 0)
    {
        for
        (
            std::initializer_list<label>::iterator iter = indexes.begin();
            iter < indexes.end();
            iter++
        )
        {
            label argIndex = std::distance(indexes.begin(), iter);
            mapIndex += (*iter)*pow(10, argSize - argIndex - 1);
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

// ************************************************************************* //
