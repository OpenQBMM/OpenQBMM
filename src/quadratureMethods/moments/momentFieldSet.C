/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2015-03-10 Alberto Passalacqua: Templated class on the type of moment and of
                                quadrature node.
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
    const autoPtr<PtrList<nodeType>>& nodes
)
:
    PtrList<momentType>
    (
        dict.lookup("moments"),
        typename momentType::iNew(distributionName, mesh, nodes)
    ),
    name_(IOobject::groupName("moments", distributionName)),
    nodes_(nodes),
    nDimensions_((*this)[0].nDimensions()),
    nMoments_((*this).size()),
    momentMap_(nMoments_)
{
    // Check on the number of moments and nodes may go here.
    // However, to keep the implementation generic, it is moved to a
    // higher-level class where the specific quadrature method is implemented.

    // Populate the moment set
    forAll(*this, mI)
    {
        momentMap_.insert
        (
            moment<momentType, nodeType>::listToLabel
            (
                this->operator[](mI).cmptOrders()
            ),
            mI
        );
    }
}


template <class momentType, class nodeType>
Foam::momentFieldSet<momentType, nodeType>::momentFieldSet
(
    const word& distributionName,
    const label nMoments,
    const autoPtr<PtrList<nodeType>>& nodes,
    const label nDimensions,
    const Map<label>& momentMap
)
:
    PtrList<momentType>(nMoments),
    name_(IOobject::groupName("moments", distributionName)),
    nodes_(nodes),
    nDimensions_(nDimensions),
    nMoments_(nMoments),
    momentMap_(momentMap)
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


// ************************************************************************* //
