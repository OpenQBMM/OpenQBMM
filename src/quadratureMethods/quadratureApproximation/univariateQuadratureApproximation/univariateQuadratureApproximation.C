/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

#include "univariateQuadratureApproximation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateQuadratureApproximation::univariateQuadratureApproximation
(
    const word& name,
    const fvMesh& mesh,
    const word support
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
    nodes_(),
    moments_(name_, *this, mesh_, nodes_),
    nPrimaryNodes_(0),
    nSecondaryNodes_(0),
    nodesNei_(),
    nodesOwn_(),
    nDimensions_(1),
    nMoments_(moments_.size()),
    momentsNei_
    (
        name_, nMoments_, nodesNei_, nDimensions_, moments_.momentMap()
    ),
    momentsOwn_
    (
        name_, nMoments_, nodesOwn_, nDimensions_, moments_.momentMap()
    ),
    momentInverter_(),
    support_(support)
{
    // Allocating nodes
    nodes_ = autoPtr<PtrList<extendedVolScalarNode> >
    (
        new PtrList<extendedVolScalarNode>
        (
            lookup("nodes"),
            Foam::extendedVolScalarNode::iNew
            (
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions(),
                moments_[0].boundaryField().types()
            )
        )
    );

    nPrimaryNodes_ = nodes_().size();
    nSecondaryNodes_ = nodes_()[0].nSecondaryNodes();

    if (nMoments_ != 2*nPrimaryNodes_ + 1)
    {
        FatalErrorInFunction
            << "Number of moments from dictionary different from number" << endl
            << "of moments calculated from primary quadrature nodes."
            << abort(FatalError);
    }

    nodesNei_ = autoPtr<PtrList<extendedSurfaceScalarNode> >
    (
        new PtrList<extendedSurfaceScalarNode>(nPrimaryNodes_)
    );

    nodesOwn_ = autoPtr<PtrList<extendedSurfaceScalarNode> >
    (
        new PtrList<extendedSurfaceScalarNode>(nPrimaryNodes_)
    );

    PtrList<extendedVolScalarNode>& nodes = nodes_();
    PtrList<extendedSurfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<extendedSurfaceScalarNode>& nodesOwn = nodesOwn_();

    // Populating interpolated nodes.
    forAll(nodes, pNodei)
    {
        extendedVolScalarNode& node(nodes[pNodei]);

        nodesNei.set
        (
            pNodei,
            new extendedSurfaceScalarNode
            (
                node.name() + "Nei",
                name_,
                nSecondaryNodes_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );

        nodesOwn.set
        (
            pNodei,
            new extendedSurfaceScalarNode
            (
                node.name() + "Own",
                name_,
                nSecondaryNodes_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );

        for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
        {

            // Commented because units of the weight would be considered twice
            // in calculations due to the product with the primary weight
            //
            //    node.secondaryWeights()[sNodei].dimensions().reset
            //    (
            //        moments_[0].dimensions();
            //    );

            node.secondaryAbscissae()[sNodei].dimensions().reset
            (
                moments_[1].dimensions()/moments_[0].dimensions()
            );
        }
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

    momentInverter_ = autoPtr<Foam::extendedMomentInversion>
    (
        Foam::extendedMomentInversion::New
        (
            subDict("extendedMomentInversionCoeff"),
            nMoments_,
            nSecondaryNodes_
        )
    );

    updateQuadrature();
    interpolateNodes();
    updateBoundaryQuadrature();
    momentsNei_.update();
    momentsOwn_.update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateQuadratureApproximation::~univariateQuadratureApproximation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::univariateQuadratureApproximation::interpolateNodes()
{
    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("nei", dimless, -1.0)
    );

    surfaceScalarField own
    (
        IOobject
        (
            "own",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("own", dimless, 1.0)
    );

    const PtrList<extendedVolScalarNode>& nodes = nodes_();
    PtrList<extendedSurfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<extendedSurfaceScalarNode>& nodesOwn = nodesOwn_();

    forAll(nodes, pNodei)
    {
        const extendedVolScalarNode& node(nodes[pNodei]);
        extendedSurfaceScalarNode& nodeOwn(nodesOwn[pNodei]);
        extendedSurfaceScalarNode& nodeNei(nodesNei[pNodei]);

        nodeOwn.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), own, "reconstruct(weight)");

        nodeOwn.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                own,
                "reconstruct(abscissa)"
            );

        nodeOwn.sigma() =
            fvc::interpolate
            (
                node.sigma(),
                own,
                "reconstruct(sigma)"
            );

        nodeNei.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), nei, "reconstruct(weight)");

        nodeNei.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                nei,
                "reconstruct(abscissa)"
            );

        nodeNei.sigma() =
            fvc::interpolate
            (
                node.sigma(),
                nei,
                "reconstruct(sigma)"
            );

        for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
        {
            // Setting interpolated secondary nodes
            nodeOwn.secondaryWeights()[sNodei] =
                fvc::interpolate
                (
                    node.secondaryWeights()[sNodei],
                    own,
                    "reconstruct(weight)"
                );

            nodeOwn.secondaryAbscissae()[sNodei] =
                fvc::interpolate
                (
                    node.secondaryAbscissae()[sNodei],
                    own,
                    "reconstruct(abscissa)"
                );

            nodeNei.secondaryWeights()[sNodei] =
                fvc::interpolate
                (
                    node.secondaryWeights()[sNodei],
                    nei,
                    "reconstruct(weight)"
                );

            nodeNei.secondaryAbscissae()[sNodei] =
                fvc::interpolate
                (
                    node.secondaryAbscissae()[sNodei],
                    nei,
                    "reconstruct(abscissa)"
                );

        }
    }
}

void Foam::univariateQuadratureApproximation::updateBoundaryQuadrature()
{
    // Recover reference to boundaryField of zero-order moment.
    // All moments will share the same BC types at a given boundary.
    const volScalarField::Boundary& bf
        = moments_().boundaryFieldRef();

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        if (m0Patch.fixesValue())
        {
            forAll(m0Patch, faceI)
            {
                univariateMomentSet momentsToInvert(nMoments_, 0, support_);

                // Copying moments from a face
                forAll(momentsToInvert, momenti)
                {
                    momentsToInvert[momenti]
                        = moments_[momenti].boundaryField()[patchi][faceI];
                }

                // Inverting them
                momentInverter_->invert(momentsToInvert);

                // Copying quadrature data to boundary face
                for (label pNodei = 0; pNodei < nPrimaryNodes_; pNodei++)
                {
                    extendedVolScalarNode& node = nodes_()[pNodei];

                    node.primaryWeight().boundaryFieldRef()[patchi][faceI]
                        = momentInverter_->primaryWeights()[pNodei];

                    node.primaryAbscissa().boundaryFieldRef()[patchi][faceI]
                        = momentInverter_->primaryAbscissae()[pNodei];

                    node.sigma().boundaryFieldRef()[patchi][faceI]
                        = momentInverter_->sigma();

                    for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
                    {
                        node.secondaryWeights()[sNodei].boundaryFieldRef()[patchi][faceI]
                            = momentInverter_->secondaryWeights()[pNodei][sNodei];

                        node.secondaryAbscissae()[sNodei].boundaryFieldRef()[patchi][faceI]
                            = momentInverter_->secondaryAbscissae()[pNodei][sNodei];
                    }
                }

            }
        }
    }
}

void Foam::univariateQuadratureApproximation::updateQuadrature()
{
    const volScalarField& m0(moments_[0]);

    PtrList<extendedVolScalarNode>& nodes(nodes_());

    forAll(m0, celli)
    {
        univariateMomentSet momentsToInvert(nMoments_, 0.0, "Gauss", support_);

        // Copying moment set from a cell to univariateMomentSet
        forAll(momentsToInvert, momenti)
        {
            momentsToInvert[momenti] = moments_[momenti][celli];
        }

        // Inverting moments and updating secondary quadrature
        momentInverter_->invert(momentsToInvert);

        // Recovering primary weights and abscissae from moment inverter
        const scalarDiagonalMatrix& pWeights(momentInverter_->primaryWeights());

        const scalarDiagonalMatrix& pAbscissae
        (
            momentInverter_->primaryAbscissae()
        );

        // Copying to fields
        for (label pNodei = 0; pNodei < nPrimaryNodes_; pNodei++)
        {
            extendedVolScalarNode& node(nodes[pNodei]);

            // Copy primary node
            node.primaryWeight()[celli] = pWeights[pNodei];
            node.primaryAbscissa()[celli] = pAbscissae[pNodei];

            // Copy secondary nodes
            PtrList<volScalarField>& sWeightFields(node.secondaryWeights());
            PtrList<volScalarField>& sAbscissaFields(node.secondaryAbscissae());

            const scalarRectangularMatrix& sWeights
            (
                momentInverter_->secondaryWeights()
            );

            const scalarRectangularMatrix& sAbscissae
            (
                momentInverter_->secondaryAbscissae()
            );

            for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
            {
                sWeightFields[sNodei][celli] = sWeights[pNodei][sNodei];
                sAbscissaFields[sNodei][celli] = sAbscissae[pNodei][sNodei];
            }

            // Copy sigma
            node.sigma()[celli] = momentInverter_->sigma();
        }
    }

    // Updating boundary conditions
    forAll(nodes, pNodei)
    {
        extendedVolScalarNode& pNode(nodes[pNodei]);

        pNode.primaryWeight().correctBoundaryConditions();
        pNode.primaryAbscissa().correctBoundaryConditions();
        pNode.sigma().correctBoundaryConditions();

        for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
        {
            pNode.secondaryWeights()[sNodei].correctBoundaryConditions();
            pNode.secondaryAbscissae()[sNodei].correctBoundaryConditions();
        }
    }

    updateBoundaryQuadrature();
    updateMoments();
}


void Foam::univariateQuadratureApproximation::updateMoments()
{
    moments_.update();
}


// ************************************************************************* //
