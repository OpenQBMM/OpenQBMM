/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    Foam::quadratureApproximation

Description
    Provides storage for the moments and the corresponding quadrature
    approximation of a distribution function.
    Methods to update the moments and the quadrature approximation are
    provided.

SourceFiles
    quadratureApproximation.C
    quadratureApproximationI.H

\*---------------------------------------------------------------------------*/

#ifndef quadratureApproximation_H
#define quadratureApproximation_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "mappedPtrList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fieldMomentInversion.H"
#include "volVelocityNode.H"
#include "surfaceVelocityNode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class univariateQuadratureApproximation Declaration
\*---------------------------------------------------------------------------*/

template <class momentType, class nodeType>
class quadratureApproximation
:
    public IOdictionary
{
protected:
    // Protected data

        typedef momentFieldSet<momentType, nodeType> momentFieldSetType;

        //- Name of the quadrature approximation
        const word name_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        const dictionary& dict_;

        //- Moment orders used to map moments
        const labelListList momentOrders_;

        //- Node indexes used to map nodes
        const labelListList nodeIndexes_;

        //- Number of nodes in each direction
        labelList nNodes_;

        //- AutoPtr to the mappedPtrList of quadrature nodes
        autoPtr<mappedPtrList<nodeType>> nodes_;

        //- List of moments
        momentFieldSetType moments_;

        //- Dimensionality of the distribution function
        const label nDimensions_;

        //- Number of moments
        const label nMoments_;

        //- Number of secondary nodes
        const label nSecondaryNodes_;

        //- Support of the distribution function
        word support_;

        //- Extended moment inversion method
        autoPtr<fieldMomentInversion> momentFieldInverter_;


public:

    //- Default name of the phase properties dictionary
    static const word propertiesName;

    // Constructors

        //- Construct from mesh
        quadratureApproximation
        (
            const word& name,
            const fvMesh& mesh,
            const word& support
        );

        //- Construct from dictionary name, name, momentFieldSet
        //  Moment fields in the momentFieldSet passed as reference are
        //  used to initialize the local copy of moments. Inversion on
        //  creation is optional.
        quadratureApproximation
        (
            const word& dictName,
            const word& name,
            const momentFieldSetType& mFieldSet,
            bool calcQuadratureOnCreation = false
        );

    //- Destructor
    virtual ~quadratureApproximation();


    // Member Functions

        // Access

            //- Return name
            inline const word name() const;

            //- Return the nmber of nodes in each direction
            const labelList& nNodes() const;

            //- Extended moment inversion method
            inline fieldMomentInversion& momentFieldInverter();

            //- Const access to the nodes
            inline const mappedPtrList<nodeType>& nodes() const;

            //- Non-Const access to the nodes
            inline mappedPtrList<nodeType>& nodes();

            //- Const access to the moments
            inline const momentFieldSetType& moments() const;

            //- Non-const access to the moments
            inline momentFieldSetType& moments();

            //- Return the number of dimensions in the distribution
            inline label nDimensions() const;

            //- Return the number of moments
            inline label nMoments() const;

            //- Return moment orders
            inline const labelListList& momentOrders() const;

            //- return node indexes
            inline const labelListList& nodeIndexes() const;


        // Edit

            //- Recalculate the quadrature nodes from the moments
            void updateQuadrature();

            //- Recalculate the quadrature nodes from the moments
            //  on boundaries
            void updateBoundaryQuadrature();
            //- Recalculate the moments from the quadrature nodes
            void updateMoments();

            //- Recalculate moments in the specified cell
            void updateLocalMoments(label celli);

            //- Recalculate the quadrature nodes from the moments in the
            //  specified cell
            bool updateLocalQuadrature
            (
                label celli,
                bool fatalErrorOnFailedRealizabilityTest = true
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "quadratureApproximationI.H"
#   include "quadratureApproximation.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
