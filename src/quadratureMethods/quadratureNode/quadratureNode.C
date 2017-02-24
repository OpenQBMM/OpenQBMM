/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2015-02-19 Alberto Passalacqua: Templated class on type of weight and abscissa.
2015-11-02 Alberto Passalacqua: Generalized initialization of fields based on
                                value_type of the field.
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

#include "quadratureNode.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template <class weightType, class abscissaType>
Foam::quadratureNode<weightType, abscissaType>::quadratureNode
(
    const word& name,
    const word& distributionName,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const dimensionSet& abscissaDimensions
)
:
    name_(IOobject::groupName(name, distributionName)),
    weight_
    (
        IOobject
        (
            IOobject::groupName(name_, "weight"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<typename weightType::value_type>
        (
            "zeroWeight",
            weightDimensions,
            pTraits<typename weightType::value_type>::zero
        )
    ),
    abscissa_
    (
        IOobject
        (
            IOobject::groupName(name_, "abscissa"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<typename abscissaType::value_type>
        (
            "zeroAbscissa",
            abscissaDimensions,
            pTraits<typename abscissaType::value_type>::zero
        )
    ),
    extended_(false)
{}


template <class weightType, class abscissaType>
Foam::quadratureNode<weightType, abscissaType>::quadratureNode
(
    const word& name,
    const word& distributionName,
    const dictionary& nodeDict,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const dimensionSet& abscissaDimensions,
    const wordList& boundaryTypes
)
:
    name_(IOobject::groupName(name, distributionName)),
    nodeDict_(nodeDict),
    weight_
    (
        IOobject
        (
            IOobject::groupName(name_, "weight"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<typename weightType::value_type>
        (
            "zeroWeight",
            weightDimensions,
            pTraits<typename weightType::value_type>::zero
        ),
        boundaryTypes
    ),
    abscissa_
    (
        IOobject
        (
            IOobject::groupName(name_, "abscissa"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<typename abscissaType::value_type>
        (
            "zeroAbscissa",
            abscissaDimensions,
            pTraits<typename abscissaType::value_type>::zero
        ),
        boundaryTypes
    ),
    extended_(false)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class weightType, class abscissaType>
Foam::quadratureNode<weightType, abscissaType>::~quadratureNode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class weightType, class abscissaType>
Foam::autoPtr<Foam::quadratureNode<weightType, abscissaType>>
Foam::quadratureNode<weightType, abscissaType>::clone() const
{
    notImplemented("quadratureNode::clone() const");
    return autoPtr<quadratureNode<weightType, abscissaType>>(NULL);
}

// ************************************************************************* //
