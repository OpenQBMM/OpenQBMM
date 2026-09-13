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
    Copyright (C) 2019-2026 Alberto Passalacqua
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
    Foam::gaussRadauMomentInversion

Description
    Performs Gauss-Radau quadrature on a moment set, which fixes a node at the
    lower end of the support.

    That abscissa is read from the dictionary of the quadrature:

    \verbatim
        univariateMomentInversion   GaussRadau;
        minKnownAbscissa            0;
    \endverbatim

    and defaults to zero, the lower end of both R+ and [0, 1]. It is written
    in the units of the abscissa.

SourceFiles
    gaussRadauMomentInversion.C
    gaussRadauMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef gaussRadauMomentInversion_H
#define gaussRadauMomentInversion_H

#include "univariateMomentInversion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class gaussRadauMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class gaussRadauMomentInversion
:
    public univariateMomentInversion
{
    // Private data

        //- Abscissa the quadrature fixes a node at
        scalar minKnownAbscissa_;

        // Bool to check if Gauss quadrature must be performed
        bool forceGauss_;


public:

    //- Runtime type information
    TypeName("GaussRadau");


    // Constructors

        //- Construct from univariateMomentSet
        gaussRadauMomentInversion
        (
            const dictionary& dict,
            const label nMaxNodes = 0
        );



    // Public member functions

        //- Modify the recurrence relation for fixed-point quadrature
        virtual void correctRecurrence
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta
        );

        //- Calculates the number of quadrature nodes
        virtual void calcNQuadratureNodes
        (
            univariateMomentSet& moments
        );

        //- Return the number of additional quadrature nodes. Gauss-Radau
        //  fixes one abscissa, which costs one node over the Gauss
        //  quadrature built from the same moments
        virtual label nAdditionalQuadraturePoints(const label nMoments) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
