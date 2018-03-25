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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class momentFieldSetType, class nodeType>
const word
Foam::quadratureApproximation<momentFieldSetType, nodeType>::name() const
{
    return name_;
}

template <class momentFieldSetType, class nodeType>
Foam::fieldMomentInversion&
Foam::quadratureApproximation<momentFieldSetType, nodeType>
::momentFieldInverter()
{
    return momentFieldInverter_();
}

template <class momentFieldSetType, class nodeType>
const Foam::mappedPtrList<nodeType>&
Foam::quadratureApproximation<momentFieldSetType, nodeType>::nodes() const
{
    return nodes_();
}

template <class momentFieldSetType, class nodeType>
Foam::mappedPtrList<nodeType>&
Foam::quadratureApproximation<momentFieldSetType, nodeType>::nodes()
{
    return nodes_();
}

template <class momentFieldSetType, class nodeType>
const momentFieldSetType&
Foam::quadratureApproximation<momentFieldSetType, nodeType>::moments() const
{
    return moments_;
}

template <class momentFieldSetType, class nodeType>
momentFieldSetType&
Foam::quadratureApproximation<momentFieldSetType, nodeType>::moments()
{
    return moments_;
}

template <class momentFieldSetType, class nodeType>
Foam::label
Foam::quadratureApproximation<momentFieldSetType, nodeType>::nDimensions() const
{
    return nDimensions_;
}

template <class momentFieldSetType, class nodeType>
Foam::label
Foam::quadratureApproximation<momentFieldSetType, nodeType>::nMoments() const
{
    return nMoments_;
}

template <class momentFieldSetType, class nodeType>
const Foam::labelListList&
Foam::quadratureApproximation<momentFieldSetType, nodeType>::momentOrders() const
{
    return momentOrders_;
}

template <class momentFieldSetType, class nodeType>
const Foam::labelListList&
Foam::quadratureApproximation<momentFieldSetType, nodeType>::nodeIndexes() const
{
    return nodeIndexes_;
}


// ************************************************************************* //
