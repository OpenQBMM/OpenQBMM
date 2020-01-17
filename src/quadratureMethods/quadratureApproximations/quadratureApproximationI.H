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
    Copyright (C) 2019 Alberto Passalacqua
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

template <class momentType, class nodeType>
const word
Foam::quadratureApproximation<momentType, nodeType>::name() const
{
    return name_;
}

template <class momentType, class nodeType>
const labelList&
Foam::quadratureApproximation<momentType, nodeType>::nNodes() const
{
    return nNodes_;
}

template <class momentType, class nodeType>
Foam::fieldMomentInversion&
Foam::quadratureApproximation<momentType, nodeType>
::momentFieldInverter()
{
    return momentFieldInverter_();
}

template <class momentType, class nodeType>
const Foam::mappedPtrList<nodeType>&
Foam::quadratureApproximation<momentType, nodeType>::nodes() const
{
    return nodes_();
}

template <class momentType, class nodeType>
Foam::mappedPtrList<nodeType>&
Foam::quadratureApproximation<momentType, nodeType>::nodes()
{
    return nodes_();
}

template <class momentType, class nodeType>
const momentFieldSet<momentType, nodeType>&
Foam::quadratureApproximation<momentType, nodeType>::moments() const
{
    return moments_;
}

template <class momentType, class nodeType>
momentFieldSet<momentType, nodeType>&
Foam::quadratureApproximation<momentType, nodeType>::moments()
{
    return moments_;
}

template <class momentType, class nodeType>
Foam::label
Foam::quadratureApproximation<momentType, nodeType>::nDimensions() const
{
    return nDimensions_;
}

template <class momentType, class nodeType>
Foam::label
Foam::quadratureApproximation<momentType, nodeType>::nMoments() const
{
    return nMoments_;
}

template <class momentType, class nodeType>
const Foam::labelListList&
Foam::quadratureApproximation<momentType, nodeType>::momentOrders() const
{
    return momentOrders_;
}

template <class momentType, class nodeType>
const Foam::labelListList&
Foam::quadratureApproximation<momentType, nodeType>::nodeIndexes() const
{
    return nodeIndexes_;
}


// ************************************************************************* //
