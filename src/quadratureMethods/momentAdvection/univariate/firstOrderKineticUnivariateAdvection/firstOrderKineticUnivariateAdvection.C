/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
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

#include "firstOrderKineticUnivariateAdvection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace univariateAdvection
{
    defineTypeNameAndDebug(firstOrderKinetic, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentAdvection,
        firstOrderKinetic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateAdvection::firstOrderKinetic::firstOrderKinetic
(
    const dictionary& dict,
    const scalarQuadratureApproximation& quadrature,
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
        name_, nMoments_, nodesNei_, nDimensions_, moments_.map(), support
    ),
    momentsOwn_
    (
        name_, nMoments_, nodesOwn_, nDimensions_, moments_.map(), support
    ),
    momentFieldInverter_()
{
    if (nMoments_ % 2 == 0)
    {
        nNodes_ = nMoments_/2;
    }
    else
    {
        nNodes_ = (nMoments_ - 1)/2 + 1;
    }

    const Map<label> map = quadrature.nodes().map();

    nodes_ = autoPtr<mappedPtrList<volScalarNode>>
    (
        new mappedPtrList<volScalarNode>(nNodes_, map)
    );

    nodesNei_ = autoPtr<mappedPtrList<surfaceScalarNode>>
    (
        new mappedPtrList<surfaceScalarNode>(nNodes_, map)
    );

    nodesOwn_ = autoPtr<mappedPtrList<surfaceScalarNode>>
    (
        new mappedPtrList<surfaceScalarNode>(nNodes_, map)
    );

    mappedPtrList<volScalarNode>& nodes = nodes_();
    mappedPtrList<surfaceScalarNode>& nodesNei = nodesNei_();
    mappedPtrList<surfaceScalarNode>& nodesOwn = nodesOwn_();

    PtrList<dimensionSet> abscissaDimensions(1);

    abscissaDimensions.set
    (
        0,
        new dimensionSet(moments_(1).dimensions()/moments_(0).dimensions())
    );

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
                moments_(0).mesh(),
                moments_(0).dimensions(),
                abscissaDimensions
            )
        );

        nodesNei.set
        (
            nodei,
            new surfaceScalarNode
            (
                "nodeRadau" + Foam::name(nodei) + "Nei",
                name_,
                moments_(0).mesh(),
                moments_(0).dimensions(),
                abscissaDimensions
            )
        );

        nodesOwn.set
        (
            nodei,
            new surfaceScalarNode
            (
                "nodeRadau" + Foam::name(nodei) + "Own",
                name_,
                moments_(0).mesh(),
                moments_(0).dimensions(),
                abscissaDimensions
            )
        );
    }

    // Setting face values of moments
    forAll(momentsNei_, momenti)
    {
        momentsNei_.set
        (
            momenti,
            new surfaceScalarMoment
            (
                name_,
                moments_(momenti).cmptOrders(),
                nodesNei_,
                fvc::interpolate(moments_(momenti)),
                "Nei"
            )
        );

        momentsOwn_.set
        (
            momenti,
            new surfaceScalarMoment
            (
                name_,
                moments_(momenti).cmptOrders(),
                nodesOwn_,
                fvc::interpolate(moments_(momenti)),
                "Own"
            )
        );
    }

    momentFieldInverter_.set
    (
        new basicFieldMomentInversion
        (
            quadrature.subDict("momentAdvection"),
            moments_[0].mesh(),
            quadrature.momentOrders(),
            quadrature.nodeIndexes(),
            nodes_()[0].velocityIndexes(),
            0
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateAdvection::firstOrderKinetic::~firstOrderKinetic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::univariateAdvection::firstOrderKinetic::interpolateNodes()
{
    const PtrList<volScalarNode>& nodes = nodes_();
    PtrList<surfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<surfaceScalarNode>& nodesOwn = nodesOwn_();

    IStringStream weightOwnLimiter("upwind");
    IStringStream abscissaOwnLimiter("upwind");

    tmp<surfaceInterpolationScheme<scalar>> weightOwnScheme
    (
        fvc::scheme<scalar>(own_, weightOwnLimiter)
    );

    tmp<surfaceInterpolationScheme<scalar>> abscissaOwnScheme
    (
        fvc::scheme<scalar>(own_, abscissaOwnLimiter)
    );

    IStringStream weightNeiLimiter("upwind");
    IStringStream abscissaNeiLimiter("upwind");

    tmp<surfaceInterpolationScheme<scalar>> weightNeiScheme
    (
        fvc::scheme<scalar>(nei_, weightNeiLimiter)
    );

    tmp<surfaceInterpolationScheme<scalar>> abscissaNeiScheme
    (
        fvc::scheme<scalar>(nei_, abscissaNeiLimiter)
    );

    forAll(nodes, rNodei)
    {
        const volScalarNode& node(nodes[rNodei]);
        surfaceScalarNode& nodeNei(nodesNei[rNodei]);
        surfaceScalarNode& nodeOwn(nodesOwn[rNodei]);

        nodeOwn.primaryWeight() =
            weightOwnScheme().interpolate(node.primaryWeight());

        nodeNei.primaryWeight() =
            weightNeiScheme().interpolate(node.primaryWeight());

        forAll(node.primaryAbscissae(), cmpt)
        {
            nodeOwn.primaryAbscissae()[cmpt] =
                abscissaOwnScheme().interpolate
                (
                    node.primaryAbscissae()[cmpt]
                );
                
            nodeNei.primaryAbscissae()[cmpt] =
                abscissaNeiScheme().interpolate
                (
                    node.primaryAbscissae()[cmpt]
                );
        }
    }
}

Foam::scalar
Foam::univariateAdvection::firstOrderKinetic::realizableCo() const
{
    // Returning 1 because the restriction of this scheme is the same CFL
    // condition of the main scheme.
    return 1.0;
}

void Foam::univariateAdvection::firstOrderKinetic::update()
{
    momentFieldInverter_().invert(moments_, nodes_());
    interpolateNodes();
    momentsNei_.update();
    momentsOwn_.update();

    dimensionedScalar zeroPhi("zero", phi_.dimensions(), Zero);

    forAll(divMoments_, divi)
    {
        divMoments_(divi) =
            fvc::surfaceIntegrate
            (
                momentsNei_[divi]*min(phi_, zeroPhi)
              + momentsOwn_[divi]*max(phi_, zeroPhi)
            );
    }
}

// ************************************************************************* //
