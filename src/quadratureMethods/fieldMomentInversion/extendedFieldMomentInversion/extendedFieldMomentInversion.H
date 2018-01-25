/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
     \\/     M anipulation  |
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
    Foam::extendedFieldMomentInversion

Description
    Inverts a volUnivariateMomentFieldSet using the extended quadrature method
    of moments.

SourceFiles
    extendedFieldMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef extendedFieldMomentInversion_H
#define extendedFieldMomentInversion_H

#include "fieldMomentInversion.H"
#include "univariateMomentSet.H"
#include "extendedMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class extendedFieldMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class extendedFieldMomentInversion
:
    public fieldMomentInversion
{
    // Private Member Functions

        //- Moment inverter
        autoPtr<extendedMomentInversion> momentInverter_;

        //- Disallow default bitwise copy construct
        extendedFieldMomentInversion(const extendedFieldMomentInversion&);

        //- Disallow default bitwise assignment
        void operator=(const extendedFieldMomentInversion&);


public:

    //- Runtime type information
    TypeName("extendedFieldMomentInversion");


    // Constructors

        //- Construct from dictionary
        extendedFieldMomentInversion
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const label nSecondaryNodes
        );


    //- Destructor
        ~extendedFieldMomentInversion();


    // Public member Functions

        // Univariate moment inversions

            //- Invert moments in the entire computational domain
            virtual void invert
            (
                const volUnivariateMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes
            );

            //- Invert moments on boundaries
            virtual void invertBoundaryMoments
            (
                const volUnivariateMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes
            );

            //- Invert moments in a single cell
            virtual bool invertLocalMoments
            (
                const volUnivariateMomentFieldSet& moments,
                mappedPtrList<volScalarNode>& nodes,
                const label celli,
                const bool fatalErrorOnFailedRealizabilityTest = true
            );


        // Multivriate moment inversions

            //- Invert moments in the entire computational domain
            virtual void invert
            (
                const volVectorMomentFieldSet& moments,
                mappedPtrList<volVectorNode>& nodes
            );

            //- Invert moments on boundaries
            virtual void invertBoundaryMoments
            (
                const volVectorMomentFieldSet& moments,
                mappedPtrList<volVectorNode>& nodes
            );

            //- Invert moments in a single cell
            virtual bool invertLocalMoments
            (
                const volVectorMomentFieldSet& moments,
                mappedPtrList<volVectorNode>& nodes,
                const label celli,
                const bool fatalErrorOnFailedRealizabilityTest = true
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
