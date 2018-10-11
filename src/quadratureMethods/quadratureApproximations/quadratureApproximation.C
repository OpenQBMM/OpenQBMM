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

#include "quadratureApproximation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class momentFieldSetType, class nodeType>
const Foam::word Foam::quadratureApproximation<momentFieldSetType, nodeType>::
propertiesName("quadratureProperties");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class momentFieldSetType, class nodeType>
Foam::quadratureApproximation<momentFieldSetType, nodeType>::
quadratureApproximation
(
    const word& name,
    const fvMesh& mesh,
    const word& support
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("quadratureProperties", name),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    mesh_(mesh),
    dict_(*this),
    momentOrders_
    (
        const_cast
        <
            const quadratureApproximation<momentFieldSetType, nodeType>&
        >(*this).lookup("moments")
    ),
    nodeIndexes_
    (
        const_cast
        <
            const quadratureApproximation<momentFieldSetType, nodeType>&
        >(*this).lookup("nodes")
    ),
    nodes_(),
    moments_(name_, *this, mesh_, nodes_, support),
    nDimensions_(moments_[0].cmptOrders().size()),
    nMoments_(moments_.size()),
    nSecondaryNodes_
    (
        lookupOrDefault<label>("nSecondaryNodes", nMoments_ + 1)
    ),
    support_(support),
    momentFieldInverter_
    (
        fieldMomentInversion::New
        (
            (*this),
            mesh_,
            momentOrders_,
            nodeIndexes_,
            nSecondaryNodes_
        )
    )
{
    // Allocating nodes
    nodes_ = autoPtr<mappedPtrList<nodeType>>
    (
        new mappedPtrList<nodeType>
        (
            lookup("nodes"),
            typename nodeType::iNew
            (
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                moments_[0].boundaryField().types(),
                momentFieldInverter_().extended(),
                nSecondaryNodes_
            )
        )
    );

    nodes_().setMap(mappedPtrList<scalar>(nodes_().size(), nodeIndexes_).map());

    updateQuadrature();
}


template<class momentFieldSetType, class nodeType>
Foam::quadratureApproximation<momentFieldSetType, nodeType>::
quadratureApproximation
(
    const word& dictName,
    const word& name,
    const momentFieldSetType& mFieldSet,
    bool calcQuadratureOnCreation
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            mFieldSet[0].mesh().time().constant(),
            mFieldSet[0].mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    name_(name),
    mesh_(mFieldSet[0].mesh()),
    dict_(*this),
    momentOrders_
    (
        const_cast
        <
            const quadratureApproximation<momentFieldSetType, nodeType>&
        >(*this).lookup("moments")
    ),
    nodeIndexes_
    (
        const_cast
        <
            const quadratureApproximation<momentFieldSetType, nodeType>&
        >(*this).lookup("nodes")
    ),
    nodes_(),
    moments_
    (
        name_,
        mFieldSet.size(),
        nodes_,
        mFieldSet.nDimensions(),
        mFieldSet.map(),
        mFieldSet.support()
    ),
    nDimensions_(mFieldSet.nDimensions()),
    nMoments_(mFieldSet.size()),
    nSecondaryNodes_
    (
        lookupOrDefault<label>("nSecondaryNodes", nMoments_ + 1)
    ),
    support_(mFieldSet.support()),
    momentFieldInverter_
    (
        fieldMomentInversion::New
        (
            (*this),
            mesh_,
            momentOrders_,
            nodeIndexes_,
            nSecondaryNodes_
        )
    )
{
    if (nSecondaryNodes_ != 0 && !momentFieldInverter_().extended())
    {
        WarningInFunction
            << "The number of secondary nodes in the quadrature" << nl
            << "    approximation is not zero, but the selected" << nl
            << "    inversion algorithm is not of extended type." << nl
            << "    Proceeding with nSecondaryNodes = 0." << nl
            << "    No extended quadrature will be computed." << nl;
    }

    forAll(moments_, mi)
    {
        moments_.set
        (
            mi,
            new volUnivariateMoment
            (
                name_ + Foam::name(mi),
                mFieldSet[mi].cmptOrders(),
                nodes_,
                mFieldSet[mi]
            )
        );
    }

    // Allocating nodes
    nodes_ = autoPtr<mappedPtrList<nodeType>>
    (
        new mappedPtrList<nodeType>
        (
            lookup("nodes"),
            typename nodeType::iNew
            (
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                moments_[0].boundaryField().types(),
                momentFieldInverter_().extended(),
                nSecondaryNodes_
            )
        )
    );

    nodes_().setMap(mappedPtrList<scalar>(nodes_().size(), nodeIndexes_).map());

    if (calcQuadratureOnCreation)
    {
        updateQuadrature();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class momentFieldSetType, class nodeType>
Foam::quadratureApproximation<momentFieldSetType, nodeType>
::~quadratureApproximation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class momentFieldSetType, class nodeType>
void Foam::quadratureApproximation<momentFieldSetType, nodeType>
::updateQuadrature()
{
    momentFieldInverter_().invert(moments_, nodes_());
    updateMoments();
}

template<class momentFieldSetType, class nodeType>
void Foam::quadratureApproximation<momentFieldSetType, nodeType>
::updateMoments()
{
    moments_.update();
}

template<class momentFieldSetType, class nodeType>
void Foam::quadratureApproximation<momentFieldSetType, nodeType>
::updateLocalMoments(label celli)
{
    moments_.updateLocalMoments(celli);
}

template<class momentFieldSetType, class nodeType>
bool Foam::quadratureApproximation<momentFieldSetType, nodeType>
::updateLocalQuadrature(label celli, bool fatalErrorOnFailedRealizabilityTest)
{
    bool realizable = momentFieldInverter_().invertLocalMoments
    (
        moments_, nodes_(), celli, false
    );

    if (!realizable && fatalErrorOnFailedRealizabilityTest)
    {
        return realizable;
    }

    moments_.updateLocalMoments(celli);

    return realizable;
}

// ************************************************************************* //
