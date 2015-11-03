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

template <class weightType, class abscissaType, class sigmaType> 
const Foam::word& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::name() const
{
    return name_;
}

template <class weightType, class abscissaType, class sigmaType>
inline bool 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::extended() const
{
    return extended_;
}

template <class weightType, class abscissaType, class sigmaType> 
Foam::label 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::nSecondaryNodes() const
{
    return nSecondaryNodes_;
}

template <class weightType, class abscissaType, class sigmaType> 
const weightType& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::primaryWeight() const
{
    return primaryWeight_;
}

template <class weightType, class abscissaType, class sigmaType> 
weightType& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::primaryWeight()
{
    return primaryWeight_;
}

template <class weightType, class abscissaType, class sigmaType> 
const abscissaType& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::primaryAbscissa() const
{
    return primaryAbscissa_;
}

template <class weightType, class abscissaType, class sigmaType> 
abscissaType& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::primaryAbscissa()
{
    return primaryAbscissa_;
}

template <class weightType, class abscissaType, class sigmaType> 
const Foam::PtrList<weightType>& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::secondaryWeights() const
{
    return secondaryWeights_;
}

template <class weightType, class abscissaType, class sigmaType> 
Foam::PtrList<weightType>& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::secondaryWeights()
{
    return secondaryWeights_;
}

template <class weightType, class abscissaType, class sigmaType> 
const Foam::PtrList<abscissaType>& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::secondaryAbscissae() const
{
    return secondaryAbscissae_;
}

template <class weightType, class abscissaType, class sigmaType> 
Foam::PtrList<abscissaType>& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::secondaryAbscissae()
{
    return secondaryAbscissae_;
}

template <class weightType, class abscissaType, class sigmaType> 
const sigmaType& 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::sigma() const
{
    return sigma_;
}

template <class weightType, class abscissaType, class sigmaType> 
sigmaType& Foam::quadratureNode<weightType, abscissaType, sigmaType>::sigma()
{
    return sigma_;
}


// ************************************************************************* //
