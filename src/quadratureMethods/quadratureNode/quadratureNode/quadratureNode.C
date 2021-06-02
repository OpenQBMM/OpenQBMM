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
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    lengthBased_(false),
    massBased_(false),
    useVolumeFraction_(false),
    secondaryWeights_(),
    secondaryAbscissae_(),
    sigmas_(),
    nSecondaryNodes_(nSecondaryNodes),
    extended_(extended)
{
    if (weightDimensions == dimless)
    {
        useVolumeFraction_ = true;
    }

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

                if (abscissaeDimensions[dimi] == dimLength)
                {
                    lengthBased_ = true;
                }
                else if (abscissaeDimensions[dimi] == dimMass)
                {
                    massBased_ = true;
                    word rhoName = IOobject::groupName("thermo:rho", name_);

                    if (mesh.foundObject<volScalarField>(rhoName))
                    {
                        rhoPtr_ = &mesh.lookupObject<volScalarField>(rhoName);
                    }
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
                    Zero
                )
            )
        );

        if (extended_)
        {
            // Allocating secondary quadrature only if the node is of extended 
            // type
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
                        dimensionedScalar("zeroWeight", dimless, Zero)
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
                            Zero
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
                    dimensionedScalar("zeroSigma", dimless, Zero)
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
    lengthBased_(false),
    massBased_(false),
    useVolumeFraction_(false),
    secondaryWeights_(),
    secondaryAbscissae_(),
    sigmas_(),
    nSecondaryNodes_(nSecondaryNodes),
    extended_(extended)
{
    if (weightDimensions == dimless)
    {
        useVolumeFraction_ = true;
    }

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

                if (abscissaeDimensions[dimi] == dimLength)
                {
                    lengthBased_ = true;
                }
                else if (abscissaeDimensions[dimi] == dimMass)
                {
                    massBased_ = true;

                    word rhoName = 
                        IOobject::groupName
                        (
                            "thermo:rho", 
                            IOobject::group(name_)
                        );

                    if (mesh.foundObject<volScalarField>(rhoName))
                    {
                        rhoPtr_ = &mesh.lookupObject<volScalarField>(rhoName);
                    }
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
                    Zero
                ),
                boundaryTypes
            )
        );

        if (extended_)
        {
            // Allocating secondary quadrature only if the node is of extended 
            // type
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
                        dimensionedScalar("zeroWeight", dimless, Zero),
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
                            Zero
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
                    dimensionedScalar("zeroSigma", dimless, Zero),
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
    return nullptr;
}


template<class scalarType, class vectorType>
Foam::tmp<Foam::volScalarField>
Foam::quadratureNode<scalarType, vectorType>::d
(
    const volScalarField& abscissa
) const
{
    if (sizeIndex_ == -1)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "d",
                    weight_.time().timeName(),
                    weight_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                weight_.mesh(),
                dimensionedScalar("d", dimLength, Zero)
            )
        );
    }

    if (lengthBased_)
    {
        return abscissa;
    }

    if (massBased_ && rhoPtr_)
    {
        return cbrt(abscissa*6.0/(Foam::constant::mathematical::pi*(*rhoPtr_)));
    }

    return cbrt(abscissa*6.0/Foam::constant::mathematical::pi);
}


template<class scalarType, class vectorType>
Foam::scalar Foam::quadratureNode<scalarType, vectorType>::d
(
    const label celli,
    const scalar& abscissa
) const
{
    if (sizeIndex_ == -1)
    {
        return 0.0;
    }

    if (lengthBased_)
    {
        return abscissa;
    }

    if (massBased_ && rhoPtr_)
    {
        return cbrt
        (
            abscissa*6.0/(Foam::constant::mathematical::pi*(*rhoPtr_)[celli])
        );
    }

    return cbrt(abscissa*6.0/Foam::constant::mathematical::pi);
}


template<class scalarType, class vectorType>
Foam::tmp<Foam::volScalarField>
Foam::quadratureNode<scalarType, vectorType>::n
(
    const volScalarField& weight,
    const volScalarField& abscissa
) const
{
    if (!useVolumeFraction_)
    {
        return weight;

    }

    tmp<volScalarField> v;

    if (massBased_ && rhoPtr_)
    {
        v = abscissa/(*rhoPtr_);
    }
    else if (lengthBased_)
    {
        v = pow3(abscissa);
    }
    else
    {
        v = abscissa;
    }

    v.ref().max(pow3(SMALL));

    return weight/v;
}


template<class scalarType, class vectorType>
Foam::scalar Foam::quadratureNode<scalarType, vectorType>::n
(
    const label celli,
    const scalar& weight,
    const scalar& abscissa
) const
{
    if (!useVolumeFraction_)
    {
        return weight;
    }

    scalar v = pow3(SMALL);

    if (massBased_ && rhoPtr_)
    {
        v = max(v, abscissa/(*rhoPtr_)[celli]);
    }
    else if (lengthBased_)
    {
        v = max(v, pow3(abscissa));
    }
    else
    {
        v = max(v, abscissa);
    }

    return weight/v;
}


// ************************************************************************* //
