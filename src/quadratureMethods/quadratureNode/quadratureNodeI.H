/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
-------------------------------------------------------------------------------
2015-02-19 Alberto Passalacqua: Templated class on type of weight and abscissa.
2015-11-05 Alberto Passalacqua: Added interface functions for secondary
                                quadrature (not implemented).
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

template <class weightType, class abscissaType>
const Foam::word& Foam::quadratureNode<weightType, abscissaType>::name() const
{
    return name_;
}

template <class weightType, class abscissaType>
inline bool Foam::quadratureNode<weightType, abscissaType>::extended() const
{
    return extended_;
}

template <class weightType, class abscissaType>
Foam::label Foam::quadratureNode<weightType, abscissaType>::
nSecondaryNodes() const
{
    NotImplemented;
    return 0;
}

template <class weightType, class abscissaType>
const weightType&
Foam::quadratureNode<weightType, abscissaType>::primaryWeight() const
{
    return weight_;
}

template <class weightType, class abscissaType>
weightType& Foam::quadratureNode<weightType, abscissaType>::primaryWeight()
{
    return weight_;
}

template <class weightType, class abscissaType>
const abscissaType&
Foam::quadratureNode<weightType, abscissaType>::primaryAbscissa() const
{
    return abscissa_;
}

template <class weightType, class abscissaType>
abscissaType& Foam::quadratureNode<weightType, abscissaType>::primaryAbscissa()
{
    return abscissa_;
}

template <class weightType, class abscissaType>
const Foam::PtrList<weightType>&
Foam::quadratureNode<weightType, abscissaType>::secondaryWeights() const
{
    NotImplemented;
    return NullObjectRef<Foam::PtrList<weightType>>();
}

template <class weightType, class abscissaType>
Foam::PtrList<weightType>&
Foam::quadratureNode<weightType, abscissaType>::secondaryWeights()
{
    NotImplemented;
    return NullObjectRef<Foam::PtrList<weightType>>();
}

template <class weightType, class abscissaType>
const Foam::PtrList<abscissaType>&
Foam::quadratureNode<weightType, abscissaType>::secondaryAbscissae() const
{
    NotImplemented;
    return NullObjectRef<Foam::PtrList<abscissaType>>();
}

template <class weightType, class abscissaType>
Foam::PtrList<abscissaType>&
Foam::quadratureNode<weightType, abscissaType>::secondaryAbscissae()
{
    NotImplemented;
    return NullObjectRef<Foam::PtrList<abscissaType>>();
}

// ************************************************************************* //
