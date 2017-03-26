/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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
    // Private member functions

        //- Disallow default bitwise copy construct
        fieldMomentInversion(const fieldMomentInversion&);

        //- Disallow default bitwise assignment
        void operator=(const fieldMomentInversion&);


protected:

    // Protected data

        //- Boolean flag to indicate extended type of inversion
        bool extended_;


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
                const label nMoments,
                const label nSecondaryNodes
            ),
            (dict, nMoments, nSecondaryNodes)
        );


    // Constructors

        //- Construct from univariateMomentSet
        fieldMomentInversion
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    //- Destructor
    virtual ~fieldMomentInversion();


    // Selectors

        static autoPtr<fieldMomentInversion> New
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    // Member Functions

        //- Invert moments in the entire computational domain
        virtual void invert
        (
            const volUnivariateMomentFieldSet& moments,
            PtrList<volScalarNode>& nodes
        ) = 0;

        //- Invert moments on boundaries
        virtual void invertBoundaryMoments
        (
            const volUnivariateMomentFieldSet& moments,
            PtrList<volScalarNode>& nodes
        ) = 0;

        //- Invert moments in a single cell
        virtual bool invertLocalMoments
        (
            const volUnivariateMomentFieldSet& moments,
            PtrList<volScalarNode>& nodes,
            const label celli,
            const bool fatalErrorOnFailedRealizabilityTest = true
        ) = 0;

        //- Return extended type
        inline bool extended();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldMomentInversionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
