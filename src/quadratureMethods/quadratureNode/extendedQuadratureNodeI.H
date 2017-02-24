/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::name() const
{
    return name_;
}

template <class weightType, class abscissaType, class sigmaType>
inline bool Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::
extended() const
{
    return extended_;
}

template <class weightType, class abscissaType, class sigmaType>
Foam::label
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::
nSecondaryNodes() const
{
    return nSecondaryNodes_;
}

template <class weightType, class abscissaType, class sigmaType>
const Foam::PtrList<weightType>&
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::
secondaryWeights() const
{
    return secondaryWeights_;
}

template <class weightType, class abscissaType, class sigmaType>
Foam::PtrList<weightType>&
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::
secondaryWeights()
{
    return secondaryWeights_;
}

template <class weightType, class abscissaType, class sigmaType>
const Foam::PtrList<abscissaType>&
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>
::secondaryAbscissae() const
{
    return secondaryAbscissae_;
}

template <class weightType, class abscissaType, class sigmaType>
Foam::PtrList<abscissaType>&
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>
::secondaryAbscissae()
{
    return secondaryAbscissae_;
}

template <class weightType, class abscissaType, class sigmaType>
const sigmaType&
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>
::sigma() const
{
    return sigma_;
}

template <class weightType, class abscissaType, class sigmaType>
sigmaType&
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::sigma()
{
    return sigma_;
}


// ************************************************************************* //
