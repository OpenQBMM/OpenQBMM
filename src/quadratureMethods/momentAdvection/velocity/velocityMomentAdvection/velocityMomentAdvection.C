/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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

#include "velocityMomentAdvection.H"
#include "IOmanip.H"
#include "wallFvPatch.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityMomentAdvection, 0);
    defineRunTimeSelectionTable(velocityMomentAdvection, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityMomentAdvection::velocityMomentAdvection
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const word& support
)
:
    name_(quadrature.name()),
    moments_(quadrature.moments()),
    nodes_(quadrature.nodes()),
    nodesNei_(),
    nodesOwn_(),
    nMoments_(moments_.size()),
    own_
    (
        IOobject
        (
            "velocityMomentAdvection:own",
            moments_(0).mesh().time().timeName(),
            moments_(0).mesh()
        ),
        moments_(0).mesh(),
        dimensionedScalar("own", dimless, 1.0)
    ),
    nei_
    (
        IOobject
        (
            "velocityMomentAdvection:nei",
            moments_(0).mesh().time().timeName(),
            moments_(0).mesh()
        ),
        moments_(0).mesh(),
        dimensionedScalar("nei", dimless, -1.0)
    ),
    support_(support),
    momentOrders_(quadrature.momentOrders()),
    nodeIndexes_(quadrature.nodeIndexes()),
    divMoments_(nMoments_),
    boundaries_(moments_[0].boundaryField().size())
{
    forAll(divMoments_, momenti)
    {
        const labelList& momentOrder = momentOrders_[momenti];

        divMoments_.set
        (
            momenti,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "divMoment"
                       + mappedList<vector>::listToWord(momentOrder),
                        name_
                    ),
                    moments_(0).mesh().time().timeName(),
                    moments_(0).mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                moments_(0).mesh(),
                dimensionedScalar
                (
                    "zero", moments_[momenti].dimensions()/dimTime, 0
                )
            )
        );
    }

    PtrList<dimensionSet> abscissaeDimensions(momentOrders_[0].size());
    labelList zeroOrder(momentOrders_[0].size(), 0);

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
    }

    nodesNei_ = autoPtr<PtrList<surfaceVelocityNode> >
    (
        new PtrList<surfaceVelocityNode>(nodes_.size())
    );

    nodesOwn_ = autoPtr<PtrList<surfaceVelocityNode> >
    (
        new PtrList<surfaceVelocityNode>(nodes_.size())
    );

    PtrList<surfaceVelocityNode>& nodesNei = nodesNei_();
    PtrList<surfaceVelocityNode>& nodesOwn = nodesOwn_();

    // Populating nodes and interpolated nodes
    forAll(nodes_, nodei)
    {
        const labelList& nodeIndex = nodeIndexes_[nodei];
        
        nodesNei.set
        (
            nodei,
            new surfaceVelocityNode
            (
                "nodeNei" + mappedList<scalar>::listToWord(nodeIndex),
                name_,
                moments_(0).mesh(),
                moments_(0).dimensions(),
                abscissaeDimensions
            )
        );

        nodesOwn.set
        (
            nodei,
            new surfaceVelocityNode
            (
                "nodeOwn" + mappedList<scalar>::listToWord(nodeIndex),
                name_,
                moments_(0).mesh(),
                moments_(0).dimensions(),
                abscissaeDimensions
            )
        );
    }

    forAll(boundaries_, patchi)
    {
        boundaries_.set
        (
            patchi,
            fvQuadraturePatch::New
            (
                moments_[0].mesh().boundary()[patchi],
                dict.optionalSubDict
                (
                    moments_[0].mesh().boundary()[patchi].name()
                ),
                quadrature,
                nodesOwn_(),
                nodesNei_()
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityMomentAdvection::~velocityMomentAdvection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityMomentAdvection::updateBoundaryConditions()
{
    forAll(boundaries_, patchi)
    {
        boundaries_[patchi].update();
    }
}

// ************************************************************************* //
