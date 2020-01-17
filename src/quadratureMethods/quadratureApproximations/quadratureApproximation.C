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
    Copyright (C) 2019 Alberto Passalacqua
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

template<class momentType, class nodeType>
const Foam::word Foam::quadratureApproximation<momentType, nodeType>::
propertiesName("quadratureProperties");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class momentType, class nodeType>
Foam::quadratureApproximation<momentType, nodeType>::
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
            const quadratureApproximation<momentType, nodeType>&
        >(*this).lookup("moments")
    ),
    nodeIndexes_
    (
        const_cast
        <
            const quadratureApproximation<momentType, nodeType>&
        >(*this).lookup("nodes")
    ),
    nNodes_(momentOrders_[0].size(), 1),
    nodes_(),
    moments_(name_, *this, mesh_, nodes_, support),
    nDimensions_(moments_[0].cmptOrders().size()),
    nMoments_(moments_.size()),
    nSecondaryNodes_
    (
        lookupOrDefault<label>("nSecondaryNodes", nMoments_ + 1)
    ),
    support_(support),
    momentFieldInverter_()
{
    forAll(nodeIndexes_, nodei)
    {
        forAll(nNodes_, dimi)
        {
            nNodes_[dimi] = max(nNodes_[dimi], nodeIndexes_[nodei][dimi] + 1);
        }
    }

    PtrList<dimensionSet> abscissaeDimensions(momentOrders_[0].size());
    labelList zeroOrder(momentOrders_[0].size(), 0);
    labelList velocityIndexes;

    forAll(abscissaeDimensions, dimi)
    {
        labelList firstOrder(zeroOrder);
        firstOrder[dimi] = 1;

        abscissaeDimensions.set
        (
            dimi,
            new dimensionSet
            (
                moments_(firstOrder).dimensions()/moments_(0).dimensions()
            )
        );

        if (abscissaeDimensions[dimi] == dimVelocity)
        {
            velocityIndexes.append(dimi);
        }
    }

    if (velocityIndexes.size() == 0)
    {
        velocityIndexes.append(-1);
    }

    momentFieldInverter_ =
        fieldMomentInversion::New
        (
            (*this),
            mesh_,
            momentOrders_,
            nodeIndexes_,
            velocityIndexes,
            nSecondaryNodes_
        );

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
                abscissaeDimensions,
                moments_[0].boundaryField().types(),
                momentFieldInverter_().extended(),
                nSecondaryNodes_
            )
        )
    );

    nodes_().setMap(mappedPtrList<scalar>(nodes_().size(), nodeIndexes_).map());

    updateQuadrature();
}


template<class momentType, class nodeType>
Foam::quadratureApproximation<momentType, nodeType>::
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
            const quadratureApproximation<momentType, nodeType>&
        >(*this).lookup("moments")
    ),
    nodeIndexes_
    (
        const_cast
        <
            const quadratureApproximation<momentType, nodeType>&
        >(*this).lookup("nodes")
    ),
    nNodes_(momentOrders_[0].size(), 1),
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
    momentFieldInverter_()
{
    forAll(nodeIndexes_, nodei)
    {
        forAll(nNodes_, dimi)
        {
            nNodes_[dimi] = max(nNodes_[dimi], nodeIndexes_[nodei][dimi] + 1);
        }
    }

    forAll(moments_, mi)
    {
        moments_.set
        (
            mi,
            new momentType
            (
                name_ + Foam::name(mi),
                mFieldSet[mi].cmptOrders(),
                nodes_,
                mFieldSet[mi]
            )
        );
    }

    PtrList<dimensionSet> abscissaeDimensions(momentOrders_[0].size());
    labelList zeroOrder(momentOrders_[0].size(), 0);
    labelList velocityIndexes;

    forAll(abscissaeDimensions, dimi)
    {
        labelList firstOrder(zeroOrder);
        firstOrder[dimi] = 1;

        abscissaeDimensions.set
        (
            dimi,
            new dimensionSet
            (
                moments_(firstOrder).dimensions()/moments_(0).dimensions()
            )
        );
        
        if (abscissaeDimensions[dimi] == dimVelocity)
        {
            velocityIndexes.append(dimi);
        }
    }

    momentFieldInverter_ =
        fieldMomentInversion::New
        (
            (*this),
            mesh_,
            momentOrders_,
            nodeIndexes_,
            velocityIndexes,
            nSecondaryNodes_
        );

    if (nSecondaryNodes_ != 0 && !momentFieldInverter_().extended())
    {
        WarningInFunction
            << "The number of secondary nodes in the quadrature" << nl
            << "    approximation is not zero, but the selected" << nl
            << "    inversion algorithm is not of extended type." << nl
            << "    Proceeding with nSecondaryNodes = 0." << nl
            << "    No extended quadrature will be computed." << nl;
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
                abscissaeDimensions,
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

template<class momentType, class nodeType>
Foam::quadratureApproximation<momentType, nodeType>
::~quadratureApproximation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class momentType, class nodeType>
void Foam::quadratureApproximation<momentType, nodeType>
::updateQuadrature()
{
    momentFieldInverter_().invert(moments_, nodes_());
    updateMoments();
}

template<class momentType, class nodeType>
void Foam::quadratureApproximation<momentType, nodeType>
::updateBoundaryQuadrature()
{
    momentFieldInverter_().invertBoundaryMoments(moments_, nodes_());
    moments_.updateBoundaries();
}

template<class momentType, class nodeType>
void Foam::quadratureApproximation<momentType, nodeType>
::updateMoments()
{
    moments_.update();
}

template<class momentType, class nodeType>
void Foam::quadratureApproximation<momentType, nodeType>
::updateLocalMoments(label celli)
{
    moments_.updateLocalMoments(celli);
}

template<class momentType, class nodeType>
bool Foam::quadratureApproximation<momentType, nodeType>
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
