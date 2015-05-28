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
    nDimensions_(1),                 
    nMoments_(2*nPrimaryNodes_ + 1),  
    momentsToInvert_(nMoments_, 0.0),
    momentInverter_
    (
        extendedMomentInversion::New
        (
            subDict("extendedMomentInversionCoeff"), 
            momentsToInvert_
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
    
    // Resetting units of quadrature nodes to make them consistent with
    // moments. Ideally this should be done at construction, but it is
    // not possible because the moment construction depends on the nodes, 
    // and nodes are not read from file
    forAll(nodes_, pNodeI)
    {
        volScalarNode& node(nodes_[pNodeI]);
        
        node.primaryWeight().dimensions().reset(moments_[0].dimensions());
        node.primaryWeight().dimensions().reset
        (
            moments_[1].dimensions()/moments_[0].dimensions()
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateQuadratureApproximation::~univariateQuadratureApproximation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::univariateQuadratureApproximation::updateQuadrature()
{
    //NOTE: Needs profiling. Allocating a lot of fields just to replace?
    //      At the same time, it reduces the number of calls to access functions
    //      significantly. Need profiling!
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
    
        sigma[cellI] = momentInverter_->sigma();
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
}


void Foam::univariateQuadratureApproximation::updateMoments()
{
    moments_.update();
}


// ************************************************************************* //
