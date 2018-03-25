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
    Foam::hyperbolicFieldMomentInversion

Description
    Invert volVectorMomentFieldSet using basic hyperbolic conditional
    moment inversion.

SourceFiles
    hyperbolicFieldMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef hyperbolicFieldMomentInversion_H
#define hyperbolicFieldMomentInversion_H

#include "fieldMomentInversion.H"
#include "multivariateMomentSet.H"
#include "hyperbolicConditionalMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class hyperbolicFieldMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class hyperbolicFieldMomentInversion
:
    public fieldMomentInversion
{
    // Private Member Functions

        //- Moment inverter
        autoPtr<hyperbolicConditionalMomentInversion> momentInverter_;

        //- Disallow default bitwise copy construct
        hyperbolicFieldMomentInversion(const hyperbolicFieldMomentInversion&);

        //- Disallow default bitwise assignment
        void operator=(const hyperbolicFieldMomentInversion&);


public:

    //- Runtime type information
    TypeName("hyperbolicFieldMomentInversion");


    // Constructors

        //- Construct from dictionary
        hyperbolicFieldMomentInversion
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const label nSecondaryNodes
        );


    //- Destructor
        ~hyperbolicFieldMomentInversion();


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
