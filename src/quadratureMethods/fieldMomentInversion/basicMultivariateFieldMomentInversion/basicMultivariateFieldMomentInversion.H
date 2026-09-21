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
    Copyright (C) 2019-2026 Alberto Passalacqua
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::basicMultivariateFieldMomentInversion

Description
    Invert the moments of a distribution of several dimensions cell by cell
    with a multivariate moment inversion, selected and set up by the
    basicMultivariateMomentInversion dictionary of the quadrature: the
    closures of velocities (CHyQMOM, CHyQMOMPlus), of a size and velocities
    (sizeCHyQMOM) and of conditional means (conditionalMeans, monoKinetic).
    The dimensions need not be velocities.

    Selected as basicMultivariateFieldMomentInversion. The names it had,
    basicVelocityFieldMomentInversion and basicVelocityMomentInversion for
    its dictionary, are still read, with a warning.

    The cells that could not be inverted, those where a conditioned value
    was clipped to its support, and those whose moments are those of no
    distribution, are counted and reported at each inversion of the mesh.
    The last of these are inverted as far as their moments go, and the
    report says by how far beyond an end of [0, 1] the canonical moment of
    the conditioning coordinate lies: the round-off the moments carry out of
    their own equations, which transport solved for each moment on its own
    leaves where a gradient is steep, or more. A canonical moment is
    divided by one minus the one below it, so the distance it is reported
    beyond an end grows with how close to an end that one is: it says which
    moment sets are worth looking at, not how far the moments are from a
    distribution.

SourceFiles
    basicMultivariateFieldMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef basicMultivariateFieldMomentInversion_H
#define basicMultivariateFieldMomentInversion_H

#include "fieldMomentInversion.H"
#include "multivariateMomentSet.H"
#include "multivariateMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class basicMultivariateFieldMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class basicMultivariateFieldMomentInversion
:
    public fieldMomentInversion
{
    // Private data

        //- Moment set used to invert the moments of a single cell or of a
        //  single boundary face. It is allocated on first use and reused, so
        //  that the inversion does not rebuild the lists and the index map of
        //  a moment set for every cell of the mesh.
        autoPtr<multivariateMomentSet> momentsToInvert_;


    // Private Member Functions

        //- Moment inverter
        autoPtr<multivariateMomentInversion> momentInverter_;

        //- Return the moment set used to invert a single cell or face,
        //  allocating it on first use
        multivariateMomentSet& localMomentSet
        (
            const volVelocityMomentFieldSet& moments
        );


public:

    //- Runtime type information
    TypeName("basicMultivariateFieldMomentInversion");


    // Constructors

        //- Construct from dictionary
        basicMultivariateFieldMomentInversion
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


        //- Disallow default bitwise copy construct
        basicMultivariateFieldMomentInversion
        (
            const basicMultivariateFieldMomentInversion&
        ) = delete;


    //- Destructor
        ~basicMultivariateFieldMomentInversion();


    // Public member Functions

        // Scalar moment inversions

            //- Invert moments in the entire computational domain
            virtual void invert
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes
            );

            //- Invert moments on boundaries
            virtual void invertBoundaryMoments
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes
            );

            //- Invert moments in a single cell
            virtual bool invertLocalMoments
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes,
                const label celli,
                const bool fatalErrorOnFailedRealizabilityTest = true
            );


        // Velocity moment inversions

            //- Invert moments in the entire computational domain
            virtual void invert
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes
            );

            //- Invert moments on boundaries
            virtual void invertBoundaryMoments
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes
            );

            //- Invert moments in a single cell
            virtual bool invertLocalMoments
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes,
                const label celli,
                const bool fatalErrorOnFailedRealizabilityTest = true
            );

        //- Whether the inversion keeps a conditioned value within its
        //  support
        virtual void clipToSupports(const bool clip)
        {
            momentInverter_->clipToSupports(clip);
        }

        //- Return the smallest accepted value for m0
        scalar smallM0() const;

        //- Return the smallest accepted value for zeta_k
        scalar smallZeta() const;
        
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicMultivariateFieldMomentInversion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
