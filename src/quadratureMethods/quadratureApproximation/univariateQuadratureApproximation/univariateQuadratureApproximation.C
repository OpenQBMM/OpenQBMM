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
    nodes_(),
    nodesRadau_(),
    moments_(name_, *this, mesh_, nodes_),
    nPrimaryNodes_(0),
    nSecondaryNodes_(0),
    nNodesRadau_(0),
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
    nodes_ = autoPtr<PtrList<extendedVolScalarNode>>
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
            << "Number of moments from dictionary different from number" << nl
            << "    of moments calculated from primary quadrature nodes."
            << abort(FatalError);
    }

    nNodesRadau_ = nPrimaryNodes_ + 1;

    nodesRadau_ = autoPtr<PtrList<basicVolScalarNode>>
    (
        new PtrList<basicVolScalarNode>(nNodesRadau_)
    );

    nodesNei_ = autoPtr<PtrList<basicSurfaceScalarNode>>
    (
        new PtrList<basicSurfaceScalarNode>(nNodesRadau_)
    );

    nodesOwn_ = autoPtr<PtrList<basicSurfaceScalarNode>>
    (
        new PtrList<basicSurfaceScalarNode>(nNodesRadau_)
    );

    PtrList<basicVolScalarNode>& nodesRadau = nodesRadau_();
    PtrList<basicSurfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<basicSurfaceScalarNode>& nodesOwn = nodesOwn_();

    // Populating Radau and interpolated nodes
    forAll(nodesRadau, rNodei)
    {
        nodesRadau.set
        (
            rNodei,
            new basicVolScalarNode
            (
                "nodeRadau" + Foam::name(rNodei),
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );

        nodesNei.set
        (
            rNodei,
            new basicSurfaceScalarNode
            (
                "nodeRadau" + Foam::name(rNodei) + "Nei",
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );

        nodesOwn.set
        (
            rNodei,
            new basicSurfaceScalarNode
            (
                "nodeRadau" + Foam::name(rNodei) + "Own",
                name_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );
    }

    // Setting face values of moments
    forAll(momentsNei_, momenti)
    {
        momentsNei_.set
        (
            momenti,
            new Foam::basicSurfaceUnivariateMoment
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
            new Foam::basicSurfaceUnivariateMoment
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

    const PtrList<basicVolScalarNode>& nodesRadau = nodesRadau_();
    PtrList<basicSurfaceScalarNode>& nodesNei = nodesNei_();
    PtrList<basicSurfaceScalarNode>& nodesOwn = nodesOwn_();

    forAll(nodesRadau, rNodei)
    {
        const basicVolScalarNode& node(nodesRadau[rNodei]);
        basicSurfaceScalarNode& nodeNei(nodesNei[rNodei]);
        basicSurfaceScalarNode& nodeOwn(nodesOwn[rNodei]);

        nodeOwn.primaryWeight() =
            fvc::interpolate(node.primaryWeight(), own, "reconstruct(weight)");

        nodeOwn.primaryAbscissa() =
            fvc::interpolate
            (
                node.primaryAbscissa(),
                own,
                "reconstruct(abscissa)"
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
    }
}

void Foam::univariateQuadratureApproximation::updateBoundaryQuadrature()
{
    // Recover reference to boundaryField of zero-order moment.
    // All moments will share the same BC types at a given boundary.
    const volScalarField::Boundary& bf = moments_().boundaryFieldRef();

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        if (m0Patch.fixesValue())
        {
            forAll(m0Patch, facei)
            {
                univariateMomentSet momentsToInvert
                (
                    nMoments_,
                    0.0,
                    "Gauss",
                    support_
                );

                univariateMomentSet momentsToInvertRadau
                (
                    nMoments_,
                    0.0,
                    "GaussRadau",
                    support_
                );

                // Copying moments from a face
                forAll(momentsToInvert, momenti)
                {
                    momentsToInvert[momenti]
                        = moments_[momenti].boundaryField()[patchi][facei];

                    momentsToInvertRadau[momenti] = momentsToInvert[momenti];
                }

                // Inverting moments for EQMOM
                momentInverter_->invert(momentsToInvert);

                // Finding Gauss-Radau quadrature
                momentsToInvertRadau.invert();

                // Copying quadrature data to boundary face
                for (label pNodei = 0; pNodei < nPrimaryNodes_; pNodei++)
                {
                    extendedVolScalarNode& node = nodes_()[pNodei];

                    node.primaryWeight().boundaryFieldRef()[patchi][facei]
                        = momentInverter_->primaryWeights()[pNodei];

                    node.primaryAbscissa().boundaryFieldRef()[patchi][facei]
                        = momentInverter_->primaryAbscissae()[pNodei];

                    node.sigma().boundaryFieldRef()[patchi][facei]
                        = momentInverter_->sigma();

                    for (label sNodei = 0; sNodei < nSecondaryNodes_; sNodei++)
                    {
                        node.secondaryWeights()[sNodei].boundaryFieldRef()[patchi][facei]
                            = momentInverter_->secondaryWeights()[pNodei][sNodei];

                        node.secondaryAbscissae()[sNodei].boundaryFieldRef()[patchi][facei]
                            = momentInverter_->secondaryAbscissae()[pNodei][sNodei];
                    }
                }

                // Copying Radau quadrature data to boundary face
                for (label rNodei = 0; rNodei < nNodesRadau_; rNodei++)
                {
                    basicVolScalarNode& nodeRadau = nodesRadau_()[rNodei];

                    if (rNodei < momentsToInvertRadau.nNodes())
                    {
                        nodeRadau.primaryWeight().boundaryFieldRef()[patchi][facei]
                            = momentsToInvertRadau.weights()[rNodei];

                        nodeRadau.primaryAbscissa().boundaryFieldRef()[patchi][facei]
                            = momentsToInvertRadau.abscissae()[rNodei];
                    }
                    else
                    {
                        nodeRadau.primaryWeight().boundaryFieldRef()[patchi][facei]
                            = 0.0;

                        nodeRadau.primaryAbscissa().boundaryFieldRef()[patchi][facei]
                            = 0.0;
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
    PtrList<basicVolScalarNode>& nodesRadau(nodesRadau_());

    forAll(m0, celli)
    {
        univariateMomentSet momentsToInvert
        (
            nMoments_,
            0.0,
            "Gauss",
            support_
        );

        univariateMomentSet momentsToInvertRadau
        (
            nMoments_,
            0.0,
            "GaussRadau",
            support_
        );

        // Copying moment set from a cell to univariateMomentSet
        forAll(momentsToInvert, momenti)
        {
            momentsToInvert[momenti] = moments_[momenti][celli];
            momentsToInvertRadau[momenti] = momentsToInvert[momenti];
        }

        // Inverting moments and updating EQMOM
        momentInverter_->invert(momentsToInvert);

        // Finding Gauss-Radau quadrature
        momentsToInvertRadau.invert();

        // Recovering primary weights and abscissae from moment inverter
        const scalarDiagonalMatrix& pWeights(momentInverter_->primaryWeights());

        const scalarDiagonalMatrix& pAbscissae
        (
            momentInverter_->primaryAbscissae()
        );

        // Recovering Gauss-Radau quadrature
        const scalarDiagonalMatrix& rWeights(momentsToInvertRadau.weights());

        const scalarDiagonalMatrix& rAbscissae
        (
            momentsToInvertRadau.abscissae()
        );

        // Copying EQMOM quadrature to fields
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

        for (label rNodei = 0; rNodei < nNodesRadau_; rNodei++)
        {
            basicVolScalarNode& node(nodesRadau[rNodei]);

            if (rNodei < momentsToInvertRadau.nNodes())
            {
                node.primaryWeight()[celli] = rWeights[rNodei];
                node.primaryAbscissa()[celli] = rAbscissae[rNodei];
            }
            else
            {
                node.primaryWeight()[celli] = 0.0;
                node.primaryAbscissa()[celli] = 0.0;
            }
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

    forAll(nodesRadau, rNodei)
    {
        basicVolScalarNode& rNode(nodesRadau[rNodei]);

        rNode.primaryWeight().correctBoundaryConditions();
        rNode.primaryAbscissa().correctBoundaryConditions();
    }

    updateBoundaryQuadrature();
    updateMoments();
}


void Foam::univariateQuadratureApproximation::updateMoments()
{
    moments_.update();
}


// ************************************************************************* //
