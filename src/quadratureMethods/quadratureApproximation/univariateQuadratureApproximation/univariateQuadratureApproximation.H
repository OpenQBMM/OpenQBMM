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

Class
    Foam::univariateQuadratureApproximation

Description
    Provides storage for the moments and the corresponding quadrature
    approximation of a univariate distribution function.
    Methods to update the moments and the quadrature approximation are
    provided.

SourceFiles
    univariateQuadratureApproximation.C
    univariateQuadratureApproximationI.H

\*---------------------------------------------------------------------------*/

#ifndef univariateQuadratureApproximation_H
#define univariateQuadratureApproximation_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "PtrList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "univariateMomentSet.H"
#include "extendedMomentInversion.H"
#include "quadratureNodes.H"
#include "moments.H"
#include "momentFieldSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class univariateQuadratureApproximation Declaration
\*---------------------------------------------------------------------------*/

class univariateQuadratureApproximation
:
    public IOdictionary
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- AutoPtr to the PtrList of quadrature nodes
        autoPtr<PtrList<extendedVolScalarNode> > nodes_;

        //- List of moments
        momentFieldSet<volUnivariateMoment, extendedVolScalarNode> moments_;
        
        //- Number of primary quadrature nodes
        label nPrimaryNodes_;

        //- Number of secondary quadrature nodes
        label nSecondaryNodes_;
        
        //- List of interpolated nodes (neighbour)
        autoPtr<PtrList<extendedSurfaceScalarNode> > nodesNei_;
        
        //- List of interpolated nodes (owner)
        autoPtr<PtrList<extendedSurfaceScalarNode> > nodesOwn_;
        
        //- Dimensionality of the distribution function
        const label nDimensions_;

        //- Number of moments
        const label nMoments_;
        
        //- List of interpolated moments (neighbour)
        momentFieldSet<surfaceUnivariateMoment, extendedSurfaceScalarNode> momentsNei_;

        //- List of interpolated moments (owner)
        momentFieldSet<surfaceUnivariateMoment, extendedSurfaceScalarNode> momentsOwn_;

        //- Moment inversion method
        autoPtr<extendedMomentInversion> momentInverter_;


    // Private member functions

        //- Calculate the m1/m0 ratio
        //tmp<volScalarField> m1ByM0();

public:

    // Constructors

        //- Construct from mesh
        explicit univariateQuadratureApproximation(const fvMesh& mesh);


    //- Destructor
    virtual ~univariateQuadratureApproximation();


    // Member Functions

        // Access

            //- Return the nodes
            inline const PtrList<extendedVolScalarNode>& nodes() const;

            //- Const access to the moments
            inline const momentFieldSet<volUnivariateMoment, extendedVolScalarNode>& 
                moments() const;

            //- Non-const access to the moments
            inline momentFieldSet<volUnivariateMoment, extendedVolScalarNode>& 
                moments();
                
            //- Const access to the interpolated moments (neighbour)
            inline const momentFieldSet
            <
                surfaceUnivariateMoment, 
                extendedSurfaceScalarNode
            >& momentsNei() const;

            //- Non-const access to the interpolated moments (neighbour)
            inline momentFieldSet
            <
                surfaceUnivariateMoment, 
                extendedSurfaceScalarNode
            >& momentsNei();
                
            //- Const access to the interpolated moments (owner)
            inline const momentFieldSet
            <
                surfaceUnivariateMoment, 
                extendedSurfaceScalarNode
            >& momentsOwn() const;

            //- Non-const access to the interpolated moments (owner)
            inline momentFieldSet
            <
                surfaceUnivariateMoment, 
                extendedSurfaceScalarNode
            >& momentsOwn();

            //- Return the number of dimensions in the distribution
            inline label nDimensions() const;

            //- Return the number of moments
            inline label nMoments() const;


        // Edit

            //- Interpolate nodes
            void interpolateNodes();
            
            //- Update quadrature nodes on fixed-value boundaries
            void updateBoundaryQuadrature();
            
            //- Recalculate the quadrature nodes from the moments
            void updateQuadrature();

            //- Recalculate the moments from the quadrature nodes
            void updateMoments();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateQuadratureApproximationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
