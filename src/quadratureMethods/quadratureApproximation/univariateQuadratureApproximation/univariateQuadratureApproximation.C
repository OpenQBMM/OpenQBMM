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
    nodes_
    ( 
        lookup("nodes"), 
        Foam::volScalarNode::iNew
        (   
            mesh_,
            dimless, //moments_[0].dimensions(),
            dimless  // moments_[1].dimensions()/moments_[0].dimensions()
        )
    ),
    moments_(*this, nodes()),
    nPrimaryNodes_(nodes_.size()),
    nSecondaryNodes_(nodes_[0].nSecondaryNodes()),
    nodesNei_(nPrimaryNodes_),
    nodesOwn_(nPrimaryNodes_),
    nDimensions_(1),                 
    nMoments_(2*nPrimaryNodes_ + 1),
    momentsNei_(nMoments_, nodesNei_, nDimensions_, moments_.momentMap()),
    momentsOwn_(nMoments_, nodesOwn_, nDimensions_, moments_.momentMap()),
    momentsToInvert_(nMoments_, 0.0),
    momentInverter_
    (
        extendedMomentInversion::New
        (
            subDict("extendedMomentInversionCoeff"), 
            momentsToInvert_,
            nSecondaryNodes_
        )
    )
{  
    if (moments_.size() != nMoments_)
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
    
    // Resetting units of nodes and populating interpolated nodes.
    forAll(nodes_, pNodeI)
    {
        volScalarNode& node(nodes_[pNodeI]);
        
        node.primaryWeight().dimensions().reset(moments_[0].dimensions());
        node.primaryAbscissa().dimensions().reset
        (
            moments_[1].dimensions()/moments_[0].dimensions()
        );
               
        nodesNei_.set
        (
            pNodeI,
            new surfaceScalarNode
            (
                node.name() + "Nei",
                nSecondaryNodes_,
                mesh_,
                moments_[0].dimensions(),
                moments_[1].dimensions()/moments_[0].dimensions()
            )
        );
        
        nodesOwn_.set
        (
            pNodeI,
            new surfaceScalarNode
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

    forAll(nodes_, pNodeI)
    {
        volScalarNode& node(nodes_[pNodeI]);
               
        nodesOwn_[pNodeI].primaryWeight() = 
            fvc::interpolate(node.primaryWeight(), own, "reconstruct(weight)");
        
        nodesOwn_[pNodeI].primaryAbscissa() = 
            fvc::interpolate
            (
                node.primaryAbscissa(), 
                own, 
                "reconstruct(abscissa)"
            );
        
        nodesOwn_[pNodeI].sigma() = 
            fvc::interpolate
            (
                node.sigma(), 
                own, 
                "reconstruct(sigma)"
            );
            
        nodesNei_[pNodeI].primaryWeight() = 
            fvc::interpolate(node.primaryWeight(), nei, "reconstruct(weight)");
        
        nodesNei_[pNodeI].primaryAbscissa() = 
            fvc::interpolate
            (
                node.primaryAbscissa(), 
                nei, 
                "reconstruct(abscissa)"
            );
            
        nodesNei_[pNodeI].sigma() = 
            fvc::interpolate
            (
                node.sigma(), 
                nei, 
                "reconstruct(sigma)"
            );
            
        for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
        {           
            // Setting interpolated secondary nodes
            nodesOwn_[pNodeI].secondaryWeights()[sNodeI] = 
                fvc::interpolate
                (
                    node.secondaryWeights()[sNodeI], 
                    own, 
                    "reconstruct(weight)"
                );
            
            nodesOwn_[pNodeI].secondaryAbscissae()[sNodeI] =
                fvc::interpolate
                (
                    node.secondaryAbscissae()[sNodeI], 
                    own, 
                    "reconstruct(abscissa)"
                );
            
            nodesNei_[pNodeI].secondaryWeights()[sNodeI] =
                fvc::interpolate
                (
                    node.secondaryWeights()[sNodeI], 
                    nei, 
                    "reconstruct(weight)"
                );
            
            nodesNei_[pNodeI].secondaryAbscissae()[sNodeI] =
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
                // Copying moments from a face
                forAll(momentsToInvert_, mI)
                {
                    momentsToInvert_[mI] 
                        = moments_[mI].boundaryField()[patchI][faceI];
                }
                
                // Inverting them
                momentInverter_->correct();
                
                // Copying quadrature data to boundary face
                for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
                {
                    volScalarNode& node = nodes_[pNodeI];
                    
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
    const label nCells = moments_().size();

    // Matrix to store primary weights
    DiagonalMatrix<scalarField> primaryWeights
    (
        nPrimaryNodes_, 
        scalarField(nCells, 0)
    );
    
    // Matrix to store primary abscissae
    DiagonalMatrix<scalarField> primaryAbscissae
    (
        nPrimaryNodes_, 
        scalarField(nCells, 0)
    );

    // Matrix to store secondary weights
    RectangularMatrix<scalarField> secondaryWeights
    (
        nPrimaryNodes_, nSecondaryNodes_, scalarField(nCells, 0)
    );
    
    // Matrix to store secondary abscissae
    RectangularMatrix<scalarField> secondaryAbscissae
    (
        nPrimaryNodes_, nSecondaryNodes_, scalarField(nCells, 0)
    );

    // Field to store sigma
    scalarField sigma(nCells, 0);

    forAll(moments_[0], cellI)
    {
        // Copying moment set from a cell to univariateMomentSet
        forAll(momentsToInvert_, mI)
        {
            momentsToInvert_[mI] = moments_[mI][cellI];
        }
       
        // Inverting moments and updating secondary quadrature
        momentInverter_->correct();

        for (label pNodeI = 0; pNodeI < nPrimaryNodes_; pNodeI++)
        {
            primaryWeights[pNodeI][cellI] 
                    = momentInverter_->primaryWeights()[pNodeI];
            
            primaryAbscissae[pNodeI][cellI] 
                    = momentInverter_->primaryAbscissae()[pNodeI];

            for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
            {
                secondaryWeights[pNodeI][sNodeI][cellI] 
                    = momentInverter_->secondaryWeights()[pNodeI][sNodeI];
            
                secondaryAbscissae[pNodeI][sNodeI][cellI] 
                        = momentInverter_->secondaryAbscissae()[pNodeI][sNodeI];
            }
        }        
    }

    forAll(nodes_, pNodeI)
    {
        // Copying primary nodes
        nodes_[pNodeI].primaryWeight().internalField().replace
        (
            0,
            primaryWeights[pNodeI]
        );

        nodes_[pNodeI].primaryWeight().correctBoundaryConditions();

        nodes_[pNodeI].primaryAbscissa().internalField().replace
        (
            0,
            primaryAbscissae[pNodeI]
        );

        nodes_[pNodeI].primaryAbscissa().correctBoundaryConditions();

        // Copying sigma
        nodes_[pNodeI].sigma().internalField().replace(0, sigma);
        nodes_[pNodeI].sigma().correctBoundaryConditions();

        // Copying secondary weights
        for (label sNodeI = 0; sNodeI < nSecondaryNodes_; sNodeI++)
        {
            nodes_[pNodeI].secondaryWeights()[sNodeI].internalField().replace
            (
                0,
                secondaryWeights[pNodeI][sNodeI]
            );

            nodes_[pNodeI].secondaryWeights()[sNodeI].correctBoundaryConditions();

            nodes_[pNodeI].secondaryAbscissae()[sNodeI].internalField().replace
            (
                0,
                secondaryAbscissae[pNodeI][sNodeI]
            );

            nodes_[pNodeI].secondaryAbscissae()[sNodeI].correctBoundaryConditions();
        }
    }
    
    updateBoundaryQuadrature();
}


void Foam::univariateQuadratureApproximation::updateMoments()
{
    moments_.update();
}


// ************************************************************************* //
