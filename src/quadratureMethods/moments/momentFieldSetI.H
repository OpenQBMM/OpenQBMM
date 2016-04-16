/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2015-03-10 Alberto Passalacqua: Templated class on the type of moment and of
                                quadrature node.
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
// template <class momentType, class nodeType>
// Foam::label
// Foam::momentFieldSet<momentType, nodeType>::dimToCmptOrder(label dim) const
// {
//     dim = nDimensions_ - dim - 1;
//
//     return pow(10, dim);
// }

template <class momentType, class nodeType>
momentType& Foam::momentFieldSet<momentType, nodeType>::operator()()
{
    return this->operator[](momentMap_[0]);
}

template <class momentType, class nodeType>
const momentType&
Foam::momentFieldSet<momentType, nodeType>::operator()() const
{
    return this->operator[](momentMap_[0]);
}

template <class momentType, class nodeType>
momentType& Foam::momentFieldSet<momentType, nodeType>::operator()(label a)
{
    return this->operator[](map(0, 0, a));
}

template <class momentType, class nodeType>
const momentType&
Foam::momentFieldSet<momentType, nodeType>::operator()(label a) const
{
    return this->operator[](map(0, 0, a));
}

template <class momentType, class nodeType>
momentType& Foam::momentFieldSet<momentType, nodeType>::operator()
(
    label a,
    label b
)
{
    return this->operator[](map(0, a, b));
}

template <class momentType, class nodeType>
const momentType& Foam::momentFieldSet<momentType, nodeType>::operator()
(
    label a,
    label b
) const
{
    return this->operator[](map(0, a, b));
}

template <class momentType, class nodeType>
momentType& Foam::momentFieldSet<momentType, nodeType>::operator()
(
    label a,
    label b,
    label c
)
{
    return this->operator[](map(a, b, c));
}

template <class momentType, class nodeType>
const momentType& Foam::momentFieldSet<momentType, nodeType>::operator()
(
    label a,
    label b,
    label c
) const
{
    return this->operator[](map(a, b, c));
}


template <class momentType, class nodeType> const Foam::Map<Foam::label>&
Foam::momentFieldSet<momentType, nodeType>::momentMap() const
{
    return momentMap_;
}


template <class momentType, class nodeType>
Foam::label Foam::momentFieldSet<momentType, nodeType>::nDimensions() const
{
    return nDimensions_;
}


// ************************************************************************* //
