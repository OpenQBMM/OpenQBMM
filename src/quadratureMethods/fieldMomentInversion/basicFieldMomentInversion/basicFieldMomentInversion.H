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
    Copyright (C) 2019 Alberto Passalacqua
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
    Foam::basicFieldMomentInversion

Description
    Invert moment a volUnivariateMomentFieldSet using basic quadrature
    inversion.

SourceFiles
    basicFieldMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef basicFieldMomentInversion_H
#define basicFieldMomentInversion_H

#include "fieldMomentInversion.H"
#include "univariateMomentSet.H"
#include "univariateMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class basicFieldMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class basicFieldMomentInversion
:
    public fieldMomentInversion
{
    // Private data

        //- Minimum value of the abscissa for fixed-point quadrature
        scalar minKnownAbscissa_;

        //- Maximum value of the abscissa for fixed-point quadrature
        scalar maxKnownAbscissa_;

        //- Number of fixed quadrature points
        label nFixedQuadraturePoints_;


    // Private Member Functions

        //- Moment inverter
        autoPtr<univariateMomentInversion> momentInverter_;


public:

    //- Runtime type information
    TypeName("basicFieldMomentInversion");


    // Constructors

        //- Construct from dictionary
        basicFieldMomentInversion
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes,
            const label nSecondaryNodes
        );

        //- Disallow default bitwise copy construct
        basicFieldMomentInversion(const basicFieldMomentInversion&) = delete;


    //- Destructor
        ~basicFieldMomentInversion();


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

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicFieldMomentInversion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
