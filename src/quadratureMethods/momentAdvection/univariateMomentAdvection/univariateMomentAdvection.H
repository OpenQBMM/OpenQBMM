/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 Alberto Passalacqua
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
    Foam::univariateMomentAdvection

Description
    Computes advection term for a moment set using a realizable advection
    scheme.

SourceFiles
    univariateMomentAdvection.C
    univariateMomentAdvectionI.H

\*---------------------------------------------------------------------------*/

#ifndef univariateMomentAdvection_H
#define univariateMomentAdvection_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "PtrList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "quadratureNodes.H"
#include "momentFieldSets.H"
#include "basicFieldMomentInversion.H"
#include "quadratureApproximations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class univariateMomentAdvection Declaration
\*---------------------------------------------------------------------------*/

class univariateMomentAdvection
{
    // Private data

        //- Name of the distribution associated to the quadrature approximation.
        const word name_;

        //- Reference to the moments to advect
        const volUnivariateMomentFieldSet& moments_;

        //- Number of moments
        const label nMoments_;

        //- Advection terms for moment transport equations
        PtrList<volScalarField> divMoments_;

        //- List of quadrature nodes used for moment advection
        autoPtr<PtrList<volScalarNode>> nodes_;

        //- List of interpolated nodes (neighbour)
        autoPtr<PtrList<surfaceScalarNode>> nodesNei_;

        //- List of interpolated nodes (owner)
        autoPtr<PtrList<surfaceScalarNode>> nodesOwn_;

        //- List of interpolated moments (neighbour)
        surfaceUnivariateMomentFieldSet momentsNei_;

        //- List of interpolated moments (owner)
        surfaceUnivariateMomentFieldSet momentsOwn_;

        //- Flux moment inversion method
        autoPtr<basicFieldMomentInversion> momentFieldInverter_;

        //- Conservative advection flux
        const surfaceScalarField phi_;

        //- Support of the distribution function
        word support_;

        //- Dimensions of the underlying distribution
        label nDimensions_;

        //- Number of nodes used for advection
        label nNodes_;

protected:

        //- Protected member functions
        void interpolateNodes();


public:

    // Constructors

        //- Construct from quadrature approximation
        univariateMomentAdvection
        (
            const univariateQuadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~univariateMomentAdvection();


    // Member Functions

        // Access

            //- Return the advection term for the moment transport equation
            inline const PtrList<volScalarField>& divMoments() const;


        // Edit

            //- Update advection
            void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentAdvectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
