/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 Alberto Passalacqua
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

#include "velocityQuadratureNode.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class scalarType, class vectorType>
Foam::velocityQuadratureNode<scalarType, vectorType>::velocityQuadratureNode
(
    const word& name,
    const word& distributionName,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const PtrList<dimensionSet>& abscissaeDimensions,
    const bool extended,
    const label nSecondaryNodes
)
:
    quadratureNode<scalarType, vectorType>
    (
        name,
        distributionName,
        mesh,
        weightDimensions,
        abscissaeDimensions,
        extended,
        nSecondaryNodes
    ),
    velocityAbscissae_(createVelocityAbscissae(this->weight_))
{}


template<class scalarType, class vectorType>
Foam::velocityQuadratureNode<scalarType, vectorType>::velocityQuadratureNode
(
    const word& name,
    const word& distributionName,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const PtrList<dimensionSet>& abscissaeDimensions,
    const wordList& boundaryTypes,
    const bool extended,
    const label nSecondaryNodes
)
:
    quadratureNode<scalarType, vectorType>
    (
        name,
        distributionName,
        mesh,
        weightDimensions,
        abscissaeDimensions,
        boundaryTypes,
        extended,
        nSecondaryNodes
    ),
    velocityAbscissae_(createVelocityAbscissae(this->weight_, boundaryTypes))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class scalarType, class vectorType>
Foam::velocityQuadratureNode<scalarType, vectorType>::
~velocityQuadratureNode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
template<class scalarType, class vectorType>
Foam::tmp<vectorType>
Foam::velocityQuadratureNode<scalarType, vectorType>::createVelocityAbscissae
(
    const scalarType& weight,
    const bool boundary
) const
{
    const fvMesh& mesh = weight.mesh();

    return tmp<vectorType>
    (
        new vectorType
        (
            IOobject
            (
                IOobject::groupName("velocityAbscissae", this->name_),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "zeroVelocityAbscissa",
                dimVelocity,
                Zero
            )
        )
    );
}*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class scalarType, class vectorType>
Foam::autoPtr<Foam::velocityQuadratureNode<scalarType, vectorType>>
Foam::velocityQuadratureNode<scalarType, vectorType>::clone() const
{
    notImplemented("velocityQuadratureNode::clone() const");
    return nullptr; //autoPtr<velocityQuadratureNode<scalarType, vectorType>>(NULL);
}

// ************************************************************************* //
