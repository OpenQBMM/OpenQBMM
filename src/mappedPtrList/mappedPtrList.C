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

#include "mappedPtrList.H"


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template <class mappedType>
Foam::word
Foam::mappedPtrList<mappedType>::listToWord(const labelList& list)
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
Foam::mappedPtrList<mappedType>::listToLabel
(
    const labelList& list,
    const label nDimensions
)
{
    label listLabel = 0;
    label size = max(nDimensions, list.size());

    forAll(list, dimi)
    {
        const label exponent = size - dimi - 1;

        // Compute 10^exponent using integer arithmetic
        label factor = 1;
        for (label i = 0; i < exponent; ++i)
        {
            factor *= 10;
        }

        listLabel += list[dimi]*factor;
    }

    return listLabel;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class mappedType> Foam::mappedPtrList<mappedType>::mappedPtrList
(
    const label size,
    const labelListList& indexes
)
:
    PtrList<mappedType>(size),
    map_(size),
    nDimensions_(0)
{
    if (indexes.size() != size)
    {
        FatalErrorInFunction
            << "Size mismatch: "
            << endl
            << "  List size = " << size
            << endl
            << "  Indexes list size = " << indexes.size()
            << exit(FatalError);
    }

    forAll(indexes, indexi)
    {
        nDimensions_ = max(nDimensions_, indexes[indexi].size());
    }

    forAll(*this, elemi)
    {
        insertKey(indexes[elemi], elemi);
    }
}


template <class mappedType> Foam::mappedPtrList<mappedType>::mappedPtrList
(
    const label size,
    const Map<label>& map
)
:
    PtrList<mappedType>(size),
    map_(map),
    nDimensions_(0)
{
    // Note: the number of dimensions is recovered from the decimal digit
    // count of the keys already in map, which undercounts leading-zero
    // components (e.g. key 1 for order {0, 1} looks one-dimensional).
    // Prefer setMap(const labelListList&) when the indexes are available.
    forAllConstIter(Map<label>, map_, iter)
    {
        label key = iter.key();
        label nD = 0;

        do
        {
            key /= 10;
            nD++;
        } while (key);

        nDimensions_ = max(nDimensions_, nD);
    }
}


template <class mappedType> Foam::mappedPtrList<mappedType>::mappedPtrList
(
    const PtrList<mappedType>& initList,
    const labelListList& indexes
)
:
    PtrList<mappedType>(initList),
    map_(initList.size()),
    nDimensions_(0)
{
    if (indexes.size() != initList.size())
    {
        FatalErrorInFunction
            << "Size mismatch: "
            << endl
            << "  List size = " << initList.size()
            << endl
            << "  Indexes list size = " << indexes.size()
            << exit(FatalError);
    }

    forAll(indexes, indexi)
    {
        nDimensions_ = max(nDimensions_, indexes[indexi].size());
    }

    forAll(*this, elemi)
    {
        insertKey(indexes[elemi], elemi);
    }
}

template <class mappedType>
template<class INew>
Foam::mappedPtrList<mappedType>::mappedPtrList(Istream& is, const INew& iNewt)
:
    PtrList<mappedType>(is, iNewt),
    nDimensions_(0)
{
    map_.resize(this->size());
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class mappedType>
Foam::label Foam::mappedPtrList<mappedType>::calcMapIndex
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
            const label argIndex = std::distance(indexes.begin(), iter);
            const label exponent = nDimensions_ - argIndex - 1;

            // Compute 10^exponent using integer arithmetic
            label factor = 1;
            for (label i = 0; i < exponent; ++i)
            {
                factor *= 10;
            }

            mapIndex += (*iter) * factor;
        }
    }

    return mapIndex;
}


template <class mappedType>
void Foam::mappedPtrList<mappedType>::insertKey
(
    const labelList& indexes,
    const label elemi
)
{
    const label key = listToLabel(indexes, nDimensions_);

    if (!map_.insert(key, elemi))
    {
        FatalErrorInFunction
            << "Duplicate mapped key " << key << " for indexes " << indexes
            << " (element " << elemi << ")." << nl
            << "Another entry already maps to this key - check for "
            << "colliding orders." << nl
            << exit(FatalError);
    }
}


template <class mappedType>
void Foam::mappedPtrList<mappedType>::setMap(const Map<label>& map)
{
    map_ = map;
    nDimensions_ = 0;

    // Note: the number of dimensions is recovered from the decimal digit
    // count of the keys already in map, which undercounts leading-zero
    // components (e.g. key 1 for order {0, 1} looks one-dimensional).
    // Prefer setMap(const labelListList&) when the indexes are available.
    forAllConstIter(Map<label>, map_, iter)
    {
        label key = iter.key();
        label nD = 0;

        do
        {
            key /= 10;
            nD++;
        } while (key);

        nDimensions_ = max(nDimensions_, nD);
    }
}


template <class mappedType>
void Foam::mappedPtrList<mappedType>::setMap(const labelListList& indexes)
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

    map_.clear();
    nDimensions_ = 0;

    forAll(indexes, indexi)
    {
        nDimensions_ = max(nDimensions_, indexes[indexi].size());
    }

    forAll(indexes, elemi)
    {
        insertKey(indexes[elemi], elemi);
    }
}


template <class mappedType>
bool Foam::mappedPtrList<mappedType>::set(const label i) const
{
    return PtrList<mappedType>::set(i);
}


template <class mappedType>
bool Foam::mappedPtrList<mappedType>::set(const labelList& list) const
{
    return PtrList<mappedType>::set(map_[listToLabel(list, nDimensions_)]);
}

template <class mappedType>
bool Foam::mappedPtrList<mappedType>::found(const labelList& list) const
{
    if (list.size() > nDimensions_)
    {
        return false;
    }

    return map_.found(listToLabel(list, nDimensions_));
}

template <class mappedType>
template <typename ...ArgsT>
bool Foam::mappedPtrList<mappedType>::found(ArgsT...args) const
{
    if (label(sizeof...(args)) > nDimensions_)
    {
        return false;
    }

    return map_.found(calcMapIndex({args...}));
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
    const labelList& list,
    mappedType* entry
)
{
    PtrList<mappedType>::set(map_[listToLabel(list, nDimensions_)], entry);
}


template <class mappedType>
void Foam::mappedPtrList<mappedType>::set
(
    const labelList& list,
    autoPtr<mappedType> entry
)
{
    PtrList<mappedType>::set
    (
        map_[listToLabel(list, nDimensions_)],
        std::move(entry)
    );
}


template <class mappedType>
void Foam::mappedPtrList<mappedType>::set
(
    const labelList& list,
    tmp<mappedType> entry
)
{
    PtrList<mappedType>::set(map_[listToLabel(list, nDimensions_)], entry);
}


template <class mappedType>
void Foam::mappedPtrList<mappedType>::setSize
(
    const label newSize,
    const labelListList& newIndexes
)
{
    if (newIndexes.size() != newSize)
    {
        FatalErrorInFunction
            << "Size mismatch: "
            << endl
            << "  New list size =" << newSize
            << endl
            << "  New indexes list size = " << newIndexes.size()
            << exit(FatalError);
    }

    // Resize the underlying list
    Foam::PtrList<mappedType>::setSize(newSize);

    // Rebuild the map with new indexes
    map_.clear();
    nDimensions_ = 0;

    forAll(newIndexes, indexi)
    {
        nDimensions_ = max(nDimensions_, newIndexes[indexi].size());
    }

    forAll(*this, elemi)
    {
        insertKey(newIndexes[elemi], elemi);
    }
}

template <class mappedType>
void Foam::mappedPtrList<mappedType>::resize
(
    const label newSize,
    const labelListList& indexes
)
{
    (*this).setSize(newSize, indexes);
}

// ************************************************************************* //
