/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
     \\/     M anipulation  |
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
template <class weightType, class abscissaType, class sigmaType>
Foam::quadratureNode<weightType, abscissaType, sigmaType>::
quadratureNode
(
    const word& name,
    const word& distributionName,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const dimensionSet& abscissaDimensions,
    const bool extended,
    const label nSecondaryNodes
)
:
    name_(IOobject::groupName(name, distributionName)),
    weight_
    (
        IOobject
        (
            IOobject::groupName("weight", name_),
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
            IOobject::groupName("abscissa", name_),
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
    secondaryWeights_(),
    secondaryAbscissae_(),
    sigma_(),
    nSecondaryNodes_(nSecondaryNodes),
    extended_(extended)
{
    if (extended_)
    {
        // Allocating secondary quadrature only if the node is of extended type
        secondaryWeights_.setSize(nSecondaryNodes_);
        secondaryAbscissae_.setSize(nSecondaryNodes_);

        // Allocating secondary weights and abscissae
        forAll(secondaryWeights_, nodei)
        {
            secondaryWeights_.set
            (
                nodei,
                new weightType
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "secondaryWeight." + Foam::name(nodei),
                            name_
                        ),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensioned<typename weightType::value_type>
                    (
                        "zeroWeight",
                        dimless,
                        pTraits<typename weightType::value_type>::zero
                    )
                )
            );

            secondaryAbscissae_.set
            (
                nodei,
                new abscissaType
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "secondaryAbscissa." + Foam::name(nodei),
                            name_
                        ),
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
                )
            );
        }

        // Allocating sigma
        sigma_ = autoPtr<sigmaType>
        (
            new sigmaType
            (
                IOobject
                (
                    IOobject::groupName("sigma", name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<typename sigmaType::value_type>
                (
                    "zeroSigma",
                    dimless,
                    pTraits<typename sigmaType::value_type>::zero
                )
            )
        );
    }
}


template <class weightType, class abscissaType, class sigmaType>
Foam::quadratureNode<weightType, abscissaType, sigmaType>::
quadratureNode
(
    const word& name,
    const word& distributionName,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const dimensionSet& abscissaDimensions,
    const wordList& boundaryTypes,
    const bool extended,
    const label nSecondaryNodes
)
:
    name_(IOobject::groupName(name, distributionName)),
    weight_
    (
        IOobject
        (
            IOobject::groupName("weight", name_),
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
            IOobject::groupName("abscissa", name_),
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
    secondaryWeights_(),
    secondaryAbscissae_(),
    sigma_(),
    nSecondaryNodes_(nSecondaryNodes),
    extended_(extended)
{
    if (extended_)
    {
        // Allocating secondary quadrature only if the node is of extended type
        secondaryWeights_.setSize(nSecondaryNodes_);
        secondaryAbscissae_.setSize(nSecondaryNodes_);

        // Allocating secondary weights and abscissae
        forAll(secondaryWeights_, nodei)
        {
            secondaryWeights_.set
            (
                nodei,
                new weightType
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "secondaryWeight." + Foam::name(nodei),
                            name_
                        ),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensioned<typename weightType::value_type>
                    (
                        "zeroWeight",
                        dimless,
                        pTraits<typename weightType::value_type>::zero
                    ),
                    boundaryTypes
                )
            );

            secondaryAbscissae_.set
            (
                nodei,
                new abscissaType
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "secondaryAbscissa." + Foam::name(nodei),
                            name_
                        ),
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
                )
            );

            sigma_ = autoPtr<sigmaType>
            (
                new sigmaType
                (
                    IOobject
                    (
                        IOobject::groupName("sigma", name_),
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensioned<typename sigmaType::value_type>
                    (
                        "zeroSigma",
                        dimless,
                        pTraits<typename sigmaType::value_type>::zero
                    ),
                    boundaryTypes
                )
            );
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType>
Foam::quadratureNode<weightType, abscissaType, sigmaType>::
~quadratureNode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType>
Foam::autoPtr
<
    Foam::quadratureNode<weightType, abscissaType, sigmaType>
>
Foam::quadratureNode<weightType, abscissaType, sigmaType>::clone() const
{
    notImplemented("quadratureNode::clone() const");
    return autoPtr
    <
        quadratureNode<weightType, abscissaType, sigmaType>
    >(NULL);
}

// ************************************************************************* //
