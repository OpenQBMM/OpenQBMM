/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "quadratureNode.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template <class weightType, class abscissaType, class sigmaType>
Foam::quadratureNode<weightType, abscissaType, sigmaType>::quadratureNode
(
    const word& name,
    const label nSecondaryNodes,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const dimensionSet& abscissaDimensions
)
:
    name_(name),
    nSecondaryNodes_(nSecondaryNodes),
    primaryWeight_
    (
        IOobject
        (
            IOobject::groupName(name_, "primaryWeight"),
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
    primaryAbscissa_
    (
        IOobject
        (
            IOobject::groupName(name_, "primaryAbscissa"),
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
    forAll(secondaryWeights_, nodeI)
    {
        secondaryWeights_.set
        (
            nodeI,
            new weightType
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        name_,
                        "secondaryWeight." + Foam::name(nodeI)
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
            nodeI,
            new abscissaType
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        name_,
                        "secondaryAbscissa." + Foam::name(nodeI)
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
Foam::quadratureNode<weightType, abscissaType, sigmaType>::quadratureNode
(
    const word& name,
    const dictionary& nodeDict,
    const fvMesh& mesh,
    const dimensionSet& weightDimensions,
    const dimensionSet& abscissaDimensions,
    const wordList& boundaryTypes
)
:
    name_(name),
    nodeDict_(nodeDict),
    nSecondaryNodes_
    (
        nodeDict_.lookupOrDefault("nSecondaryNodes", 10)
    ),
    primaryWeight_
    (
        IOobject
        (
            IOobject::groupName(name_, "primaryWeight"),
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
    primaryAbscissa_
    (
        IOobject
        (
            IOobject::groupName(name_, "primaryAbscissa"),
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
    forAll(secondaryWeights_, nodeI)
    {
        secondaryWeights_.set
        (
            nodeI,
            new weightType
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        name_,
                        "secondaryWeight." + Foam::name(nodeI)
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
            nodeI,
            new abscissaType
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        name_,
                        "secondaryAbscissa." + Foam::name(nodeI)
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
Foam::quadratureNode<weightType, abscissaType, sigmaType>::~quadratureNode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType> 
Foam::autoPtr<Foam::quadratureNode<weightType, abscissaType, sigmaType> > 
Foam::quadratureNode<weightType, abscissaType, sigmaType>::clone() const
{
    notImplemented("quadratureNode::clone() const");
    return autoPtr<quadratureNode<weightType, abscissaType, sigmaType> >(NULL);
}

// ************************************************************************* //
