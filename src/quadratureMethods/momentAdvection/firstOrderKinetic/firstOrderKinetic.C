/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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

#include "firstOrderKinetic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(firstOrderKinetic, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentAdvection,
        firstOrderKinetic,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::firstOrderKinetic::firstOrderKinetic
(
    const dictionary& dict,
    const univariateQuadratureApproximation& quadrature,
    const surfaceScalarField& phi,
    const word& support
)
:
    univariateMomentAdvection(dict, quadrature, phi, support),
    nodes_(),
    nodesNei_(),
    nodesOwn_(),
    momentsNei_
    (
        name_, nMoments_, nodesNei_, nDimensions_, moments_.momentMap(), support
    ),
    momentsOwn_
    (
        name_, nMoments_, nodesOwn_, nDimensions_, moments_.momentMap(), support
    ),
    momentFieldInverter_
    (
        new basicFieldMomentInversion
        (
            quadrature.subDict("momentAdvection"),
            nMoments_,
            0
        )
    )
{
    if (nMoments_ % 2 == 0)
    {
        nNodes_ = nMoments_/2;
    }
    else
    {
        nNodes_ = (nMoments_ - 1)/2 + 1;
    }

    nodes_ = autoPtr<PtrList<volScalarNode>>
    (
        new PtrList<volScalarNode>(nNodes_)
    );

    nodesNei_ = autoPtr<PtrList<surfaceScalarNode>>
    (
        new PtrList<surfaceScalarNode>(nNodes_)
    );

    nodesOwn_ = autoPtr<PtrList<surfaceScalarNode>>
    (
        new PtrList<surfaceScalarNode>(nNodes_)
    );

    PtrList<volScalarNode>& nodes = nodes_();
    PtrList<surfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<surfaceScalarNode>& nodesOwn = nodesOwn_();

    // Populating nodes and interpolated nodes
    forAll(nodes, nodei)
    {
        nodes.set
        (
            nodei,
            new volScalarNode
            (
                "nodeAdvection" + Foam::name(nodei),
                name_,
                moments_[0].mesh(),
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                false
            )
        );

        nodesNei.set
        (
            nodei,
            new surfaceScalarNode
            (
                "nodeRadau" + Foam::name(nodei) + "Nei",
                name_,
                moments_[0].mesh(),
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                false
            )
        );

        nodesOwn.set
        (
            nodei,
            new surfaceScalarNode
            (
                "nodeRadau" + Foam::name(nodei) + "Own",
                name_,
                moments_[0].mesh(),
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                false
            )
        );
    }

    // Setting face values of moments
    forAll(momentsNei_, momenti)
    {
        momentsNei_.set
        (
            momenti,
            new Foam::surfaceUnivariateMoment
            (
                name_,
                moments_[momenti].cmptOrders(),
                nodesNei_,
                fvc::interpolate(moments_[momenti])
            )
        );

        momentsOwn_.set
        (
            momenti,
            new Foam::surfaceUnivariateMoment
            (
                name_,
                moments_[momenti].cmptOrders(),
                nodesOwn_,
                fvc::interpolate(moments_[momenti])
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::firstOrderKinetic::~firstOrderKinetic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::firstOrderKinetic::interpolateNodes()
{
    const PtrList<volScalarNode>& nodes = nodes_();
    PtrList<surfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<surfaceScalarNode>& nodesOwn = nodesOwn_();

    forAll(nodes, rNodei)
    {
        const volScalarNode& node(nodes[rNodei]);
        surfaceScalarNode& nodeNei(nodesNei[rNodei]);
        surfaceScalarNode& nodeOwn(nodesOwn[rNodei]);

        nodeOwn.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), own_, "reconstruct(weight)");

        nodeOwn.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                own_,
                "reconstruct(abscissa)"
            );

        nodeNei.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), nei_, "reconstruct(weight)");

        nodeNei.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                nei_,
                "reconstruct(abscissa)"
            );
    }
}

Foam::scalar Foam::firstOrderKinetic::realizableCo()
{
    // Returning 1 because the restriction of this scheme is the same CFL
    // condition of the main scheme.
    return 1.0;
}

void Foam::firstOrderKinetic::update()
{
    momentFieldInverter_().invert(moments_, nodes_());
    interpolateNodes();
    momentsNei_.update();
    momentsOwn_.update();

    dimensionedScalar zeroPhi("zero", phi_.dimensions(), 0.0);

    forAll(divMoments_, divi)
    {
        volScalarField divMoment
        (
            IOobject
            (
                "divMoment",
                moments_[0].mesh().time().timeName(),
                moments_[0].mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            moments_[0].mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        );

        surfaceScalarField mFlux
        (
            momentsNei_[divi]*min(phi_, zeroPhi)
          + momentsOwn_[divi]*max(phi_, zeroPhi)
        );

        fvc::surfaceIntegrate(divMoment.ref(), mFlux);
        divMoment.ref().dimensions().reset(moments_[divi].dimensions()/dimTime);

        divMoments_[divi].replace(0, divMoment);
    }
}

// ************************************************************************* //
