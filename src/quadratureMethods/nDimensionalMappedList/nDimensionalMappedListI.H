/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class mappedType>
mappedType& Foam::nDimensionalMappedList<mappedType>::operator()()
{
    return this->operator[](map_[0]);
}

template <class mappedType>
const mappedType&
Foam::nDimensionalMappedList<mappedType>::operator()() const
{
    return this->operator[](map_[0]);
}

template <class mappedType>
mappedType& Foam::nDimensionalMappedList<mappedType>::operator()(label a)
{
    return this->operator[](map(0, 0, 0, 0, a));
}

template <class mappedType>
const mappedType&
Foam::nDimensionalMappedList<mappedType>::operator()(label a) const
{
    return this->operator[](map(0, 0, 0, 0, a));
}

template <class mappedType>
mappedType& Foam::nDimensionalMappedList<mappedType>::operator()
(
    label a,
    label b
)
{
    return this->operator[](map(0, 0, 0, a, b));
}

template <class mappedType>
const mappedType& Foam::nDimensionalMappedList<mappedType>::operator()
(
    label a,
    label b
) const
{
    return this->operator[](map(0, 0, 0, a, b));
}

template <class mappedType>
mappedType& Foam::nDimensionalMappedList<mappedType>::operator()
(
    label a,
    label b,
    label c
)
{
    return this->operator[](map(0, 0, a, b, c));
}

template <class mappedType>
const mappedType& Foam::nDimensionalMappedList<mappedType>::operator()
(
    label a,
    label b,
    label c
) const
{
    return this->operator[](map(0, 0, a, b, c));
}

template <class mappedType>
mappedType& Foam::nDimensionalMappedList<mappedType>::operator()
(
    label a,
    label b,
    label c,
    label d
)
{
    return this->operator[](map(0, a, b, c, d));
}

template <class mappedType>
const mappedType& Foam::nDimensionalMappedList<mappedType>::operator()
(
    label a,
    label b,
    label c,
    label d
) const
{
    return this->operator[](map(0, a, b, c, d));
}

template <class mappedType>
mappedType& Foam::nDimensionalMappedList<mappedType>::operator()
(
    label a,
    label b,
    label c,
    label d,
    label e
)
{
    return this->operator[](map(a, b, c, d, e));
}

template <class mappedType>
const mappedType& Foam::nDimensionalMappedList<mappedType>::operator()
(
    label a,
    label b,
    label c,
    label d,
    label e
) const
{
    return this->operator[](map(a, b, c, d, e));
}

template <class mappedType>
const mappedType&
Foam::nDimensionalMappedList<mappedType>::operator()
(
    const Foam::labelList& l
) const
{
    return this->operator[](map_[listToLabel(l)]);
}

template <class mappedType>
mappedType&
Foam::nDimensionalMappedList<mappedType>::operator()
(
    const Foam::labelList& l
)
{
    return this->operator[](map_[listToLabel(l)]);
}

template <class mappedType> const Foam::Map<Foam::label>&
Foam::nDimensionalMappedList<mappedType>::map() const
{
    return map_;
}


template <class mappedType>
Foam::label Foam::nDimensionalMappedList<mappedType>::nDims() const
{
    return nDims_;
}


// ************************************************************************* //
