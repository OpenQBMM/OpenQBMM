/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2012-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
-------------------------------------------------------------------------------
2015-03-10 Alberto Passalacqua: Templated class on the type of moment and of
                                quadrature node.
2017-03-26 Alberto Passalacqua: Added the capability to recompute moments
                                locally.
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

#include "momentFieldSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class momentType, class nodeType>
Foam::momentFieldSet<momentType, nodeType>::momentFieldSet
(
    const word& distributionName,
    const dictionary& dict,
    const fvMesh& mesh,
    const autoPtr<mappedPtrList<nodeType>>& nodes,
    const word& support
)
:
    mappedPtrList<momentType>
    (
        dict.lookup("moments"),
        typename momentType::iNew(distributionName, mesh, nodes)
    ),
    name_(IOobject::groupName("moments", distributionName)),
    nodes_(nodes),
    nDimensions_((*this)[0].nDimensions()),
    nMoments_((*this).size()),
    support_(support)
{
    Map<label> momentMap(nMoments_);

    // Populate the moment set
    forAll(*this, mI)
    {
        momentMap.insert
        (
            moment<momentType, nodeType>::listToLabel
            (
                this->operator[](mI).cmptOrders()
            ),
            mI
        );
    }
    this->setMap(momentMap);
}


template <class momentType, class nodeType>
Foam::momentFieldSet<momentType, nodeType>::momentFieldSet
(
    const word& distributionName,
    const label nMoments,
    const autoPtr<mappedPtrList<nodeType>>& nodes,
    const label nDimensions,
    const Map<label>& momentMap,
    const word& support
)
:
    mappedPtrList<momentType>(nMoments, momentMap),
    name_(IOobject::groupName("moments", distributionName)),
    nodes_(nodes),
    nDimensions_(nDimensions),
    nMoments_(nMoments),
    support_(support)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class momentType, class nodeType>
Foam::momentFieldSet<momentType, nodeType>::~momentFieldSet()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class momentType, class nodeType>
void Foam::momentFieldSet<momentType, nodeType>::update()
{
    forAll(*this, mI)
    {
        this->operator[](mI).update();
    }
}

template <class momentType, class nodeType>
void Foam::momentFieldSet<momentType, nodeType>::updateBoundaries()
{
    forAll(*this, mI)
    {
        this->operator[](mI).updateBoundaries();
    }
}

template <class momentType, class nodeType>
void Foam::momentFieldSet<momentType, nodeType>
::updateLocalMoments(label elemi)
{
    forAll(*this, mI)
    {
        this->operator[](mI).updateLocalMoment(elemi);
    }
}

// ************************************************************************* //
