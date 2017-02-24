/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

#include "extendedQuadratureNode.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template <class weightType, class abscissaType, class sigmaType>
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::
extendedQuadratureNode
(
    const word& name,
    const word& distributionName,
    const label nSecondaryNodes,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const dimensionSet& abscissaDimensions
)
:
    quadratureNode<weightType, abscissaType>
    (
        name,
        distributionName,
        mesh,
        weightDimensions,
        abscissaDimensions
    ),
    name_(IOobject::groupName(name, distributionName)),
    nSecondaryNodes_(nSecondaryNodes),
    secondaryWeights_(nSecondaryNodes_),
    secondaryAbscissae_(nSecondaryNodes_),
    sigma_
    (
        IOobject
        (
            IOobject::groupName(name_, "sigma"),
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
    ),
    extended_(true)
{
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
                        name_,
                        "secondaryWeight." + Foam::name(nodei)
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
                        name_,
                        "secondaryAbscissa." + Foam::name(nodei)
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
}


template <class weightType, class abscissaType, class sigmaType>
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::
extendedQuadratureNode
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
    quadratureNode<weightType, abscissaType>
    (
        name,
        distributionName,
        nodeDict,
        mesh,
        weightDimensions,
        abscissaDimensions,
        boundaryTypes
    ),
    name_(IOobject::groupName(name, distributionName)),
    nSecondaryNodes_
    (
        nodeDict.lookupOrDefault("nSecondaryNodes", 10)
    ),
    secondaryWeights_(nSecondaryNodes_),
    secondaryAbscissae_(nSecondaryNodes_),
    sigma_
    (
        IOobject
        (
            IOobject::groupName(name_, "sigma"),
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
    ),
    extended_(true)
{
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
                        name_,
                        "secondaryWeight." + Foam::name(nodei)
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
                        name_,
                        "secondaryAbscissa." + Foam::name(nodei)
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
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType>
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::
~extendedQuadratureNode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType>
Foam::autoPtr
<
    Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>
>
Foam::extendedQuadratureNode<weightType, abscissaType, sigmaType>::clone() const
{
    notImplemented("extendedQuadratureNode::clone() const");
    return autoPtr
    <
        extendedQuadratureNode<weightType, abscissaType, sigmaType>
    >(NULL);
}

// ************************************************************************* //
