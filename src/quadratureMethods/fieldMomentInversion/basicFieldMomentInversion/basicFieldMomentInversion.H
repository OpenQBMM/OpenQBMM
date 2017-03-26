/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

        //- Disallow default bitwise copy construct
        basicFieldMomentInversion(const basicFieldMomentInversion&);

        //- Disallow default bitwise assignment
        void operator=(const basicFieldMomentInversion&);


public:

    //- Runtime type information
    TypeName("basicFieldMomentInversion");


    // Constructors

        //- Construct from dictionary
        basicFieldMomentInversion
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    //- Destructor
        ~basicFieldMomentInversion();


    // Public member Functions

        //- Invert moments in the entire computational domain
        virtual void invert
        (
            const volUnivariateMomentFieldSet& moments,
            PtrList<volScalarNode>& nodes
        );

        //- Invert moments on boundaries
        virtual void invertBoundaryMoments
        (
            const volUnivariateMomentFieldSet& moments,
            PtrList<volScalarNode>& nodes
        );

        //- Invert moments in a single cell
        virtual bool invertLocalMoments
        (
            const volUnivariateMomentFieldSet& moments,
            PtrList<volScalarNode>& nodes,
            const label celli,
            const bool fatalErrorOnFailedRealizabilityTest = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
