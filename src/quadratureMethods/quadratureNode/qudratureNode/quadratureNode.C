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
template<class scalarType, class vectorType>
Foam::quadratureNode<scalarType, vectorType>::quadratureNode
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
    abscissae_(),
    scalarIndexes_(),
    velocityIndexes_(),
    sizeIndex_(-1),
    massBased_(false),
    secondaryWeights_(),
    secondaryAbscissae_(),
    sigmas_(),
    nSecondaryNodes_(nSecondaryNodes),
    extended_(extended)
{
    forAll(abscissaeDimensions, dimi)
    {
        if (abscissaeDimensions[dimi] == dimVelocity)
        {
            velocityIndexes_.append(dimi);
        }
        else
        {
            scalarIndexes_.append(dimi);

            if
            (
                (abscissaeDimensions[dimi] == dimMass)
             || (abscissaeDimensions[dimi] == dimVolume)
             || (abscissaeDimensions[dimi] == dimLength)
            )
            {
                if (sizeIndex_ != -1)
                {
                    FatalErrorInFunction
                        << "Only one abscissae can be sized based."
                        << abort(FatalError);
                }
                sizeIndex_ = dimi;

                if
                (
                    (abscissaeDimensions[dimi] == dimMass)
                 || (abscissaeDimensions[dimi] == dimVolume)
                )
                {
                    massBased_ = true;
                }
            }
        }
    }

    abscissae_.resize(scalarIndexes_.size());
    if (extended_)
    {
        secondaryWeights_.resize(scalarIndexes_.size());
        secondaryAbscissae_.resize(scalarIndexes_.size());
        sigmas_.resize(scalarIndexes_.size());
    }

    forAll(abscissae_, dimi)
    {
        label cmpt = scalarIndexes_[dimi];
        abscissae_.set
        (
            dimi,
            new abscissaType
            (
                IOobject
                (
                    IOobject::groupName("abscissa" + Foam::name(dimi), name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "zeroAbscissa",
                    abscissaeDimensions[cmpt],
                    0.0
                )
            )
        );

        if (extended_)
        {
            // Allocating secondary quadrature only if the node is of extended type
            secondaryWeights_.set
            (
                dimi,
                new PtrList<weightType>(nSecondaryNodes_)
            );
            secondaryAbscissae_.set
            (
                dimi,
                new abscissaeType(nSecondaryNodes_)
            );

            // Allocating secondary weights and abscissae
            forAll(secondaryWeights_[dimi], sNodei)
            {
                secondaryWeights_[dimi].set
                (
                    sNodei,
                    new weightType
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "secondaryWeight"
                              + Foam::name(dimi)
                              + '.'
                              + Foam::name(sNodei),
                                name_
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

                secondaryAbscissae_[dimi].set
                (
                    sNodei,
                    new abscissaType
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "secondaryAbscissa"
                              + Foam::name(dimi)
                              + '.'
                              + Foam::name(sNodei),
                                name_
                            ),
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "zeroAbscissa",
                            abscissaeDimensions[cmpt],
                            0.0
                        )
                    )
                );
            }

            // Allocating sigma
            sigmas_.set
            (
                dimi,
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
                    dimensionedScalar("zeroSigma", dimless, 0.0)
                )
            );
        }
    }
}


template<class scalarType, class vectorType>
Foam::quadratureNode<scalarType, vectorType>::quadratureNode
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
        ),
        boundaryTypes
    ),
    abscissae_(),
    scalarIndexes_(),
    velocityIndexes_(),
    sizeIndex_(-1),
    massBased_(false),
    secondaryWeights_(),
    secondaryAbscissae_(),
    sigmas_(),
    nSecondaryNodes_(nSecondaryNodes),
    extended_(extended)
{
    forAll(abscissaeDimensions, dimi)
    {
        if (abscissaeDimensions[dimi] == dimVelocity)
        {
            velocityIndexes_.append(dimi);
        }
        else
        {
            scalarIndexes_.append(dimi);

            if
            (
                (abscissaeDimensions[dimi] == dimMass)
             || (abscissaeDimensions[dimi] == dimVolume)
             || (abscissaeDimensions[dimi] == dimLength)
            )
            {
                if (sizeIndex_ != -1)
                {
                    FatalErrorInFunction
                        << "Only one abscissae can be sized based."
                        << abort(FatalError);
                }
                sizeIndex_ = dimi;

                if
                (
                    (abscissaeDimensions[dimi] == dimMass)
                 || (abscissaeDimensions[dimi] == dimVolume)
                )
                {
                    massBased_ = true;
                }
            }
        }
    }

    abscissae_.resize(scalarIndexes_.size());
    if (extended_)
    {
        secondaryWeights_.resize(scalarIndexes_.size());
        secondaryAbscissae_.resize(scalarIndexes_.size());
        sigmas_.resize(scalarIndexes_.size());
    }

    forAll(abscissae_, dimi)
    {
        label cmpt = scalarIndexes_[dimi];
        abscissae_.set
        (
            dimi,
            new abscissaType
            (
                IOobject
                (
                    IOobject::groupName("abscissa" + Foam::name(dimi), name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "zeroAbscissa",
                    abscissaeDimensions[cmpt],
                    0.0
                ),
                boundaryTypes
            )
        );

        if (extended_)
        {
            // Allocating secondary quadrature only if the node is of extended type
            secondaryWeights_.set
            (
                dimi,
                new PtrList<weightType>(nSecondaryNodes_)
            );
            secondaryAbscissae_.set
            (
                dimi,
                new abscissaeType(nSecondaryNodes_)
            );

            // Allocating secondary weights and abscissae
            forAll(secondaryWeights_[dimi], sNodei)
            {
                secondaryWeights_[dimi].set
                (
                    sNodei,
                    new weightType
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "secondaryWeight"
                              + Foam::name(dimi)
                              + '.'
                              + Foam::name(sNodei),
                                name_
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

                secondaryAbscissae_[dimi].set
                (
                    sNodei,
                    new abscissaType
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                "secondaryAbscissa"
                              + Foam::name(dimi)
                              + '.'
                              + Foam::name(sNodei),
                                name_
                            ),
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh,
                        dimensionedScalar
                        (
                            "zeroAbscissa",
                            abscissaeDimensions[cmpt],
                            0.0
                        ),
                    boundaryTypes
                    )
                );
            }

            // Allocating sigma
            sigmas_.set
            (
                dimi,
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
                    dimensionedScalar("zeroSigma", dimless, 0.0),
                    boundaryTypes
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class scalarType, class vectorType>
Foam::quadratureNode<scalarType, vectorType>::~quadratureNode()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class scalarType, class vectorType>
Foam::autoPtr<Foam::quadratureNode<scalarType, vectorType>>
Foam::quadratureNode<scalarType, vectorType>::clone() const
{
    notImplemented("quadratureNode::clone() const");
    return autoPtr<quadratureNode<scalarType, vectorType>>(NULL);
}

// ************************************************************************* //
