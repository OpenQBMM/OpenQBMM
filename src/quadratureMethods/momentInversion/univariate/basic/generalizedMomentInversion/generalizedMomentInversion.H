/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 Alberto Passalacqua
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
    Foam::generalizedMomentInversion

Description
    Performs generalized QMOM (GQMOM) quadrature on a set of 2N+1 moments with
    exponent nu.

SourceFiles
    generalizedMomentInversion.C
    generalizedMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef generalizedMomentInversion_H
#define generalizedMomentInversion_H

#include "univariateMomentInversion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class generalizedMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class generalizedMomentInversion
:
    public univariateMomentInversion
{
    // Private data

        // Exponent
        scalar nu_;

        // NDF type for R+ support
        word ndfTypeRPlus_;

        // Maximum number of nodes for generalized quadrature
        label nMaxNodes_;

        // Number of nodes from regular quadrature (internally calculated)
        label nRegularQuadratureNodes_;

        // Number of additional quadrature nodes for generalized quadrature
        label nAdditionalQuadratureNodes_;


    // Private member functions

        // Correct recurrence relation for R support
        void correctRecurrenceR
        (
            scalarList& alpha,
            scalarList& beta
        );

        // Correct recurrence relation for R+ support
        void correctRecurrenceRPlus
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta
        );

        // Correct recurrence relation for [0, 1] support
        void correctRecurrence01
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta
        );

public:

    //- Runtime type information
    TypeName("GQMOM");


    // Constructors

        //- Construct from univariateMomentSet
        generalizedMomentInversion
        (
            const dictionary& dict,
            const label nMaxNodes = 0
        );


    //- Destructor
    virtual ~generalizedMomentInversion();


    // Public member functions

        //- Modify the recurrence relation for fixed-point quadrature
        virtual void correctRecurrence
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta,
            const scalar minKnownAbscissa = 0,
            const scalar maxKnownAbscissa = 0
        );

        //- Calculates the number of quadrature nodes
        virtual void calcNQuadratureNodes
        (
            univariateMomentSet& moments
        );

        //- Invert the set of moments to compute weights and abscissae
        void invert
        (
            univariateMomentSet& moments,
            const scalar minKnownAbscissa = 0,
            const scalar maxKnownAbscissa = 0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
