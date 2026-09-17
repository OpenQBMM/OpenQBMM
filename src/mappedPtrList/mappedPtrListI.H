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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class mappedType>
template <typename ...ArgsT>
mappedType& Foam::mappedPtrList<mappedType>::operator()(ArgsT ...args)
{
    const label elemi = index_.find({args...});

    return
        this->operator[]
        (
            elemi < 0 ? existingPosition(labelList({args...})) : elemi
        );
}

template <class mappedType>
template <typename ...ArgsT>
const mappedType&
Foam::mappedPtrList<mappedType>::operator()(ArgsT ...args) const
{
    const label elemi = index_.find({args...});

    return
        this->operator[]
        (
            elemi < 0 ? existingPosition(labelList({args...})) : elemi
        );
}

template <class mappedType>
const mappedType& Foam::mappedPtrList<mappedType>::operator()
(
    const Foam::labelList& list
) const
{
    return this->operator[](existingPosition(list));
}

template <class mappedType>
mappedType& Foam::mappedPtrList<mappedType>::operator()
(
    const Foam::labelList& list
)
{
    return this->operator[](existingPosition(list));
}

template <class mappedType>
Foam::labelListList Foam::mappedPtrList<mappedType>::indexes() const
{
    return index_.indexes();
}

template <class mappedType>
Foam::label Foam::mappedPtrList<mappedType>::nDimensions() const noexcept
{
    return index_.nDimensions();
}

// ************************************************************************* //
