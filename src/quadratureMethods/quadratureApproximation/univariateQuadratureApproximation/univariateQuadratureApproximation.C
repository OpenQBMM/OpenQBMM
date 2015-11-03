/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "quadratureProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    nodes_(),
    moments_(*this, mesh_, nodes_),
    nPrimaryNodes_(0),
    nSecondaryNodes_(0),
    nodesNei_(),
    nodesOwn_(),
    nDimensions_(1),                 
    nMoments_(moments_.size()),
    momentsNei_(nMoments_, nodesNei_, nDimensions_, moments_.momentMap()),
    momentsOwn_(nMoments_, nodesOwn_, nDimensions_, moments_.momentMap()),
    momentInverter_()
{  
    // Allocating nodes
    nodes_ = autoPtr<PtrList<extendedVolScalarNode> >
    (
        new PtrList<extendedVolScalarNode>
        (
            lookup("nodes"), 
            Foam::extendedVolScalarNode::iNew
            (   
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
        FatalErrorIn
        (
            "Foam::univariateQuadratureApproximation::univariateQuadratureApproximation\n"
            "(\n"
            "    const fvMesh& mesh\n"
            ")"
        )   << "Number of moments from dictionary different from number\n"
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
    forAll(nodes, pNodeI)
    {
        extendedVolScalarNode& node(nodes[pNodeI]);
              
        nodesNei.set
        (
            pNodeI,
            new extendedSurfaceScalarNode
            (
                node.name() + "Nei",
                nSecondaryNodes_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );

        nodesOwn.set
        (
            pNodeI,
            new extendedSurfaceScalarNode
            (
                node.name() + "Own",
                nSecondaryNodes_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );

        for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
        {
        
            // Commented because units of the weight would be considered twice
            // in calculations due to the product with the primary weight
            //
            //    node.secondaryWeights()[sNodeI].dimensions().reset
            //    (
            //        moments_[0].dimensions();
            //    );

            node.secondaryAbscissae()[sNodeI].dimensions().reset
            (
                moments_[1].dimensions()/moments_[0].dimensions()
            );
        }
    }
    
    // Setting face values of moments
    forAll(momentsNei_, mI)
    {
        momentsNei_.set
        (
            mI,
            new Foam::surfaceUnivariateMoment
            (
                moments_[mI].cmptOrders(),
                nodesNei_,
                fvc::interpolate(moments_[mI])
            )
        );
        
        momentsOwn_.set
        (
            mI,
            new Foam::surfaceUnivariateMoment
            (
                moments_[mI].cmptOrders(),
                nodesOwn_,
                fvc::interpolate(moments_[mI])
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
    
    forAll(nodes, pNodeI)
    {
        const extendedVolScalarNode& node(nodes[pNodeI]);
        extendedSurfaceScalarNode& nodeOwn(nodesOwn[pNodeI]);
        extendedSurfaceScalarNode& nodeNei(nodesNei[pNodeI]);
               
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
            
        for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
        {           
            // Setting interpolated secondary nodes
            nodeOwn.secondaryWeights()[sNodeI] = 
                fvc::interpolate
                (
                    node.secondaryWeights()[sNodeI], 
                    own, 
                    "reconstruct(weight)"
                );
            
            nodeOwn.secondaryAbscissae()[sNodeI] =
                fvc::interpolate
                (
                    node.secondaryAbscissae()[sNodeI], 
                    own, 
                    "reconstruct(abscissa)"
                );
            
            nodeNei.secondaryWeights()[sNodeI] =
                fvc::interpolate
                (
                    node.secondaryWeights()[sNodeI], 
                    nei, 
                    "reconstruct(weight)"
                );
            
            nodeNei.secondaryAbscissae()[sNodeI] =
                fvc::interpolate
                (
                    node.secondaryAbscissae()[sNodeI], 
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
    volScalarField::GeometricBoundaryField& bf = moments_().boundaryField();
    
    forAll(bf, patchI)
    {
        fvPatchScalarField& m0Patch = bf[patchI];
        
        if (m0Patch.fixesValue())
        {
            forAll(m0Patch, faceI)
            {
                univariateMomentSet momentsToInvert(nMoments_, 0);

                // Copying moments from a face
                forAll(momentsToInvert, mI)
                {
                    momentsToInvert[mI] 
                        = moments_[mI].boundaryField()[patchI][faceI];
                }

                // Inverting them
                momentInverter_->invert(momentsToInvert);

                // Copying quadrature data to boundary face
                for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
                {
                    extendedVolScalarNode& node = nodes_()[pNodeI];
                    
                    node.primaryWeight().boundaryField()[patchI][faceI]
                        = momentInverter_->primaryWeights()[pNodeI];
            
                    node.primaryAbscissa().boundaryField()[patchI][faceI] 
                        = momentInverter_->primaryAbscissae()[pNodeI];
                        
                    node.sigma().boundaryField()[patchI][faceI]
                        = momentInverter_->sigma();

                    for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
                    {
                        node.secondaryWeights()[sNodeI].boundaryField()[patchI][faceI] 
                            = momentInverter_->secondaryWeights()[pNodeI][sNodeI];
            
                        node.secondaryAbscissae()[sNodeI].boundaryField()[patchI][faceI]
                            = momentInverter_->secondaryAbscissae()[pNodeI][sNodeI];
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

    forAll(m0, cellI)
    {
        univariateMomentSet momentsToInvert(nMoments_, 0.0);

        // Copying moment set from a cell to univariateMomentSet
        forAll(momentsToInvert, mI)
        {
            momentsToInvert[mI] = moments_[mI][cellI];
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
        for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
        {
            extendedVolScalarNode& node(nodes[pNodeI]);

            // Copy primary node
            node.primaryWeight()[cellI] = pWeights[pNodeI];
            node.primaryAbscissa()[cellI] = pAbscissae[pNodeI];
            
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

            for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
            {
                sWeightFields[sNodeI][cellI] = sWeights[pNodeI][sNodeI];
                sAbscissaFields[sNodeI][cellI] = sAbscissae[pNodeI][sNodeI];
            }

            // Copy sigma
            node.sigma()[cellI] = momentInverter_->sigma();
        }        
    }

    // Updating boundary conditions
    forAll(nodes, pNodeI)
    {
        extendedVolScalarNode& pNode(nodes[pNodeI]);
        
        pNode.primaryWeight().correctBoundaryConditions();
        pNode.primaryAbscissa().correctBoundaryConditions();
        pNode.sigma().correctBoundaryConditions();

        for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
        {
            pNode.secondaryWeights()[sNodeI].correctBoundaryConditions();
            pNode.secondaryAbscissae()[sNodeI].correctBoundaryConditions();
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
