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
    Copyright (C) 2019-2025 Alberto Passalacqua
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
    Foam::fieldMomentInversion

Description
    Abstract class for run-time selection of field moment inversion.

SourceFiles
    fieldMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef fieldMomentInversion_H
#define fieldMomentInversion_H

#include "runTimeSelectionTables.H"
#include "quadratureNodes.H"
#include "multivariateMomentSet.H"
#include "univariateMomentSet.H"
#include "momentFieldSets.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class fieldMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class fieldMomentInversion
{
protected:

    // Protected data

        //- Const refrence to the mesh
        const fvMesh& mesh_;

        //- Moment orders used to map moments
        const labelListList& momentOrders_;

        //- Node indexes used to map nodes
        const labelListList& nodeIndexes_;

        //- Whether the population of the domain has been reported as close
        //  to the threshold a cell is taken as empty at
        mutable bool warnedSmallM0_;


    // Protected member functions

        //- Warn, once, when the largest zero-order moment the field reaches
        //  is close enough to smallM0 that cells holding a real part of the
        //  population would be taken as empty
        void checkSmallM0(const volScalarField& m0) const;


public:

    //- Runtime type information
    TypeName("fieldMomentInversion");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fieldMomentInversion,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const labelListList& momentOrders,
                const labelListList& nodeIndexes,
                const labelList& velocityIndexes
            ),
            (
                dict,
                mesh,
                momentOrders,
                nodeIndexes,
                velocityIndexes
            )
        );


    // Constructors

        //- Construct from univariateMomentSet
        fieldMomentInversion
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );

        //- Disallow default bitwise copy construct
        fieldMomentInversion(const fieldMomentInversion&) = delete;


    //- Destructor
    virtual ~fieldMomentInversion();


    // Selectors

        //- The type of field inversion a quadrature dictionary selects. The
        //  one of several dimensions was named
        //  basicVelocityFieldMomentInversion when its dimensions could only
        //  be a size and velocities; the old name is read as the new one,
        //  with a warning.
        static word fieldInversionType(const dictionary& dict);

        //- The dictionary of the inversion of a quadrature of several
        //  dimensions, basicMultivariateMomentInversion, found under the
        //  name it had, basicVelocityMomentInversion, with a warning
        static const dictionary& multivariateInversionDict
        (
            const dictionary& dict
        );

        //- The orders of the moments a case carries and the indexes of
        //  its nodes, from the quadrature dictionary: the inversion it
        //  selects owns them, and the case gives the counts its closure
        //  leaves open. A case that lists its moments and nodes itself is
        //  taken as it is.
        static void momentSet
        (
            const dictionary& dict,
            labelListList& momentOrders,
            labelListList& nodeIndexes
        );

        //- The support of each dimension of the moments of a quadrature of
        //  nDimensions dimensions: for a velocity quadrature, the size, if
        //  it has one, on the positive reals and the velocity on the reals;
        //  for a univariate one, the support the model gives, which its
        //  physics fixes
        static List<supportType> momentSupports
        (
            const dictionary& dict,
            const label nDimensions,
            const List<supportType>& modelSupports
        );

        //- The moments and the nodes the inversion the dictionary selects
        //  builds, and whether it builds them at all
        static bool buildMomentSet
        (
            const dictionary& dict,
            const word& fieldInversion,
            labelListList& momentOrders,
            labelListList& nodeIndexes
        );

        static autoPtr<fieldMomentInversion> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    // Public member Functions

        // Scalar moment inversions

            //- Invert moments in the entire computational domain
            virtual void invert
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes
            ) = 0;

            //- Invert moments on boundaries
            virtual void invertBoundaryMoments
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes
            ) = 0;

            //- Invert moments in a single cell
            virtual bool invertLocalMoments
            (
                const volScalarMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes,
                const label celli,
                const bool fatalErrorOnFailedRealizabilityTest = true
            ) = 0;


        // Velocity moment inversions

            //- Invert moments in the entire computational domain
            virtual void invert
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes
            ) = 0;

            //- Invert moments on boundaries
            virtual void invertBoundaryMoments
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes
            ) = 0;

            //- Invert moments in a single cell
            virtual bool invertLocalMoments
            (
                const volVelocityMomentFieldSet& moments,
                mappedPtrList<volVelocityNode>& nodes,
                const label celli,
                const bool fatalErrorOnFailedRealizabilityTest = true
            ) = 0;

        //- Return the smallest accepted value for m0
        virtual scalar smallM0() const = 0;

        //- Return the smallest accepted value for m0
        virtual scalar smallZeta() const = 0;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldMomentInversion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
