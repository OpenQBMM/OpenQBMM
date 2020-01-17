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
    Copyright (C) 2019 Alberto Passalacqua
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
    fieldMomentInversionI.H

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

        //- Boolean flag to indicate extended type of inversion
        bool extended_;

        //- Moment orders used to map moments
        const labelListList& momentOrders_;

        //- Node indexes used to map nodes
        const labelListList& nodeIndexes_;


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
                const labelList& velocityIndexes,
                const label nSecondaryNodes
            ),
            (
                dict,
                mesh,
                momentOrders,
                nodeIndexes,
                velocityIndexes,
                nSecondaryNodes
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
            const labelList& velocityIndexes,
            const label nSecondaryNodes
        );

        //- Disallow default bitwise copy construct
        fieldMomentInversion(const fieldMomentInversion&) = delete;


    //- Destructor
    virtual ~fieldMomentInversion();


    // Selectors

        static autoPtr<fieldMomentInversion> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes,
            const label nSecondaryNodes
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

        //- Return extended type
        inline bool extended();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fieldMomentInversion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldMomentInversionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
