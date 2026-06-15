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

template<class scalarType, class vectorType>
const Foam::word&
Foam::quadratureNode<scalarType, vectorType>::name() const
{
    return name_;
}

template<class scalarType, class vectorType>
const Foam::labelList&
Foam::quadratureNode<scalarType, vectorType>::scalarIndexes() const
{
    return scalarIndexes_;
}

template<class scalarType, class vectorType>
const Foam::labelList&
Foam::quadratureNode<scalarType, vectorType>::velocityIndexes() const
{
    return velocityIndexes_;
}

template<class scalarType, class vectorType>
Foam::label
Foam::quadratureNode<scalarType, vectorType>::sizeIndex() const
{
    return sizeIndex_;
}

template<class scalarType, class vectorType>
bool Foam::quadratureNode<scalarType, vectorType>::lengthBased() const
{
    return lengthBased_;
}

template<class scalarType, class vectorType>
bool Foam::quadratureNode<scalarType, vectorType>::useVolumeFraction() const
{
    return useVolumeFraction_;
}

template<class scalarType, class vectorType>
const scalarType&
Foam::quadratureNode<scalarType, vectorType>::weight() const
{
    return weight_;
}

template<class scalarType, class vectorType>
scalarType&
Foam::quadratureNode<scalarType, vectorType>::weight()
{
    return weight_;
}

template<class scalarType, class vectorType>
const Foam::PtrList<scalarType>&
Foam::quadratureNode<scalarType, vectorType>::abscissae() const
{
    return abscissae_;
}

template<class scalarType, class vectorType>
Foam::PtrList<scalarType>&
Foam::quadratureNode<scalarType, vectorType>::abscissae()
{
    return abscissae_;
}

template<class scalarType, class vectorType>
const vectorType&
Foam::quadratureNode<scalarType, vectorType>::velocityAbscissae() const
{
    FatalErrorInFunction
        << "Attempt to access velocity abscissae of a quadrature node" << nl
        << "    with no velocity abscissae. "
        << abort(FatalError);

    return const_cast<const vectorType&>(velocityAbscissaeType::null());
}

template<class scalarType, class vectorType>
vectorType&
Foam::quadratureNode<scalarType, vectorType>::velocityAbscissae()
{
    FatalErrorInFunction
        << "Attempt to access velocity abscissae of a quadrature node" << nl
        << "    with no velocity abscissae. "
        << abort(FatalError);
    return const_cast<vectorType&>(velocityAbscissaeType::null());
}

// ************************************************************************* //
