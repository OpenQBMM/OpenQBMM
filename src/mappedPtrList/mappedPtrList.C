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

#include "mappedPtrList.H"


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template <class mappedType>
Foam::word
Foam::mappedPtrList<mappedType>::listToWord(const labelList& lst)
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
Foam::mappedPtrList<mappedType>::listToLabel(const labelList& lst)
{
    label l = 0;

    forAll(lst, dimi)
    {
        l += lst[dimi]*pow(10, lst.size() - dimi - 1);
    }

    return l;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class mappedType> Foam::mappedPtrList<mappedType>::mappedPtrList
(
    const label size,
    const labelListList& indexes
)
:
    PtrList<mappedType>(size),
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

template <class mappedType> Foam::mappedPtrList<mappedType>::mappedPtrList
(
    const label size,
    const Map<label>& map
)
:
    PtrList<mappedType>(size),
    map_(map)
{}

template <class mappedType> Foam::mappedPtrList<mappedType>::mappedPtrList
(
    const PtrList<mappedType>& initList,
    const labelListList& indexes
)
:
    PtrList<mappedType>(initList),
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

template <class mappedType>
template<class INew>
Foam::mappedPtrList<mappedType>::mappedPtrList(Istream& is, const INew& iNewt)
:
    PtrList<mappedType>(is, iNewt)
{
    map_.resize(this->size());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class mappedType>
Foam::mappedPtrList<mappedType>::~mappedPtrList()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class mappedType>
Foam::label Foam::mappedPtrList<mappedType>::calcMapIndex
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
void Foam::mappedPtrList<mappedType>::setMap(const Map<label>& map)
{
    map_ = map;
}

template <class mappedType>
void Foam::mappedPtrList<mappedType>::set
(
    const label i,
    mappedType* entry
)
{
    PtrList<mappedType>::set(i, entry);
}

template <class mappedType>
void Foam::mappedPtrList<mappedType>::set
(
    const labelList& l,
    mappedType* entry
)
{
    PtrList<mappedType>::set(map_[listToLabel(l)], entry);
}

template <class mappedType>
void Foam::mappedPtrList<mappedType>::set
(
    const labelList& l,
    autoPtr<mappedType> entry
)
{
    PtrList<mappedType>::set(map_[listToLabel(l)], entry);
}

template <class mappedType>
void Foam::mappedPtrList<mappedType>::set
(
    const labelList& l,
    tmp<mappedType> entry
)
{
    PtrList<mappedType>::set(map_[listToLabel(l)], entry);
}
// ************************************************************************* //