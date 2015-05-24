/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2015-03-09 Alberto Passalacqua: Templated class on the type of field used to
                                store the moment and on the type of quadrature
                                node.
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

template <class fieldType, class nodeType> 
const Foam::word& Foam::moment<fieldType, nodeType>::name() const
{
    return name_;
}

template <class fieldType, class nodeType> 
const Foam::word& Foam::moment<fieldType, nodeType>::keyword() const
{
    return name();
}

template <class fieldType, class nodeType> 
Foam::label Foam::moment<fieldType, nodeType>::nDimensions() const
{
    return nDimensions_;
}

template <class fieldType, class nodeType> 
Foam::label Foam::moment<fieldType, nodeType>::order() const
{
    return order_;
}

template <class fieldType, class nodeType> 
const Foam::labelList& 
Foam::moment<fieldType, nodeType>::cmptOrders() const
{
    return cmptOrders_;
}


// ************************************************************************* //
