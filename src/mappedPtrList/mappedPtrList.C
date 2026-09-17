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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class mappedType>
Foam::label Foam::mappedPtrList<mappedType>::existingPosition
(
    const labelList& list
) const
{
    const label elemi = index_.find(list);

    if (elemi < 0)
    {
        FatalErrorInFunction
            << "No element of orders " << list << " in a mapped list of "
            << index_.nDimensions() << " dimensions." << nl
            << exit(FatalError);
    }

    return elemi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class mappedType> Foam::mappedPtrList<mappedType>::mappedPtrList
(
    const label size,
    const labelListList& indexes
)
:
    PtrList<mappedType>(size),
    index_()
{
    setMap(indexes);
}


template <class mappedType> Foam::mappedPtrList<mappedType>::mappedPtrList
(
    const PtrList<mappedType>& initList,
    const labelListList& indexes
)
:
    PtrList<mappedType>(initList),
    index_()
{
    setMap(indexes);
}


template <class mappedType>
template<class INew>
Foam::mappedPtrList<mappedType>::mappedPtrList(Istream& is, const INew& iNewt)
:
    PtrList<mappedType>(is, iNewt),
    index_()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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

    index_ = mappedListKey(indexes);
}

template <class mappedType>
bool Foam::mappedPtrList<mappedType>::set(const label i) const
{
    return PtrList<mappedType>::set(i);
}

template <class mappedType>
bool Foam::mappedPtrList<mappedType>::set(const labelList& list) const
{
    return PtrList<mappedType>::set(existingPosition(list));
}

template <class mappedType>
bool Foam::mappedPtrList<mappedType>::found(const labelList& list) const
{
    // Access takes zeros past the dimensions of the list as the element
    // without them; found is asked whether the orders, as given, belong
    // to the list
    return list.size() <= index_.nDimensions() && index_.find(list) >= 0;
}

template <class mappedType>
template <typename ...ArgsT>
bool Foam::mappedPtrList<mappedType>::found(ArgsT...args) const
{
    return
        label(sizeof...(args)) <= index_.nDimensions()
     && index_.find({args...}) >= 0;
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
    PtrList<mappedType>::set(existingPosition(list), entry);
}

template <class mappedType>
void Foam::mappedPtrList<mappedType>::set
(
    const labelList& list,
    autoPtr<mappedType> entry
)
{
    PtrList<mappedType>::set(existingPosition(list), std::move(entry));
}

template <class mappedType>
void Foam::mappedPtrList<mappedType>::set
(
    const labelList& list,
    tmp<mappedType> entry
)
{
    PtrList<mappedType>::set(existingPosition(list), entry);
}

template <class mappedType>
void Foam::mappedPtrList<mappedType>::setSize
(
    const label newSize,
    const labelListList& newIndexes
)
{
    Foam::PtrList<mappedType>::setSize(newSize);
    setMap(newIndexes);
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
