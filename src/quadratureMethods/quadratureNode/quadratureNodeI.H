/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
-------------------------------------------------------------------------------
2015-02-19 Alberto Passalacqua: Templated class on type of weight and abscissa.
2015-03-08 Alberto Passalacqua: Generalised implementation to include secondary
                                quadrature and the sigma parameter of kernel 
                                density functions used in the extended 
                                quadrature method of moments.
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
Foam::label 
Foam::quadratureNode<weightType, abscissaType>::nSecondaryNodes() const
{
    return nSecondaryNodes_;
}

template <class weightType, class abscissaType> 
const weightType& 
Foam::quadratureNode<weightType, abscissaType>::primaryWeight() const
{
    return primaryWeight_;
}

template <class weightType, class abscissaType> 
weightType& Foam::quadratureNode<weightType, abscissaType>::primaryWeight()
{
    return primaryWeight_;
}

template <class weightType, class abscissaType> 
const abscissaType& 
Foam::quadratureNode<weightType, abscissaType>::primaryAbscissa() const
{
    return primaryAbscissa_;
}

template <class weightType, class abscissaType> 
abscissaType& Foam::quadratureNode<weightType, abscissaType>::primaryAbscissa()
{
    return primaryAbscissa_;
}

template <class weightType, class abscissaType> 
const Foam::PtrList<weightType>& 
Foam::quadratureNode<weightType, abscissaType>::secondaryWeights() const
{
    return secondaryWeights_;
}

template <class weightType, class abscissaType> 
Foam::PtrList<weightType>& 
Foam::quadratureNode<weightType, abscissaType>::secondaryWeights()
{
    return secondaryWeights_;
}

template <class weightType, class abscissaType> 
const Foam::PtrList<abscissaType>& 
Foam::quadratureNode<weightType, abscissaType>::secondaryAbscissae() const
{
    return secondaryAbscissae_;
}

template <class weightType, class abscissaType> 
Foam::PtrList<abscissaType>& 
Foam::quadratureNode<weightType, abscissaType>::secondaryAbscissae()
{
    return secondaryAbscissae_;
}

template <class weightType, class abscissaType> 
const weightType& Foam::quadratureNode<weightType, abscissaType>::sigma() const
{
    return sigma_;
}

template <class weightType, class abscissaType> 
weightType& Foam::quadratureNode<weightType, abscissaType>::sigma()
{
    return sigma_;
}


// ************************************************************************* //
