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
template <class weightType, class abscissaType>
Foam::quadratureNode<weightType, abscissaType>::quadratureNode
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
        dimensionedScalar("zeroWeight", weightDimensions, 0.0)
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
        dimensionedScalar("zeroAbscissa", abscissaDimensions, 0.0)
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
        dimensionedScalar("zeroSigma", dimless, 0.0)
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
                dimensionedScalar("zeroWeight", dimless, 0.0)
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
                dimensionedScalar("zeroAbscissa", abscissaDimensions, 0.0)            
            )
        );
    }
}


template <class weightType, class abscissaType>
Foam::quadratureNode<weightType, abscissaType>::quadratureNode
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
    nSecondaryNodes_(nodeDict_.lookupOrDefault("nSecondaryNodes", 7)),
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
        dimensionedScalar("zeroWeight", weightDimensions, 0.0),
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
        dimensionedScalar("zeroAbscissa", abscissaDimensions, 0.0),
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
        dimensionedScalar("zeroSigma", dimless, 0.0),
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
                dimensionedScalar("zeroWeight", dimless, 0.0),
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
                dimensionedScalar("zeroAbscissa", abscissaDimensions, 0.0),
                boundaryTypes
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class weightType, class abscissaType>
Foam::quadratureNode<weightType, abscissaType>::~quadratureNode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class weightType, class abscissaType> 
Foam::autoPtr<Foam::quadratureNode<weightType, abscissaType> > 
Foam::quadratureNode<weightType, abscissaType>::clone() const
{
    notImplemented("quadratureNode::clone() const");
    return autoPtr<quadratureNode<weightType, abscissaType> >(NULL);
}

// ************************************************************************* //
