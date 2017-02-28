/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 Alberto Passalacqua
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
    Foam::gammaEQMOM

Description
    Extended quadrature method of moments with gamma kernel density functions.

    References
    \verbatim
        "An extended quadrature method of moments for population balance
        equations"
        C Yuan, F Laurent, R O Fox
        Journal of Aerosol Science
        Volume 51, Pages 1-23, 2012
    \endverbatim

    \verbatim
        "The theory of canonical moments with applications in Statistics,
        Probability and Analysis"
        H Dette, W J Studden
        Wiley & Sons, 1997
    \endverbatim

SourceFiles
    gammaEQMOM.C

\*---------------------------------------------------------------------------*/

#ifndef gammaEQMOM_H
#define gammaEQMOM_H

#include "extendedMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class gammaEQMOM Declaration
\*---------------------------------------------------------------------------*/

class gammaEQMOM
:
    public extendedMomentInversion
{
public:

    //- Runtime type information
    TypeName("gamma");


    // Constructors

        //- Construct from a dictionary and label
        gammaEQMOM
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    //- Destructor
    virtual ~gammaEQMOM();


    // Member Functions

        //- Computes the secondary abscissa
        virtual scalar secondaryAbscissa
        (
            scalar primaryAbscissa,
            scalar secondaryAbscissa,
            scalar sigma
        );

        //- Compute moments from starred moments
        virtual void momentsStarToMoments
        (
            scalar sigma,
            univariateMomentSet& moments,
            const univariateMomentSet& momentsStar
        );

        //- Compute the starred moments
        virtual void momentsToMomentsStar
        (
            scalar sigma,
            const univariateMomentSet& moments,
            univariateMomentSet& momentsStar
        );

        //- Compute the last moment from starred moments
        virtual scalar m2N(scalar sigma, univariateMomentSet momentsStar);

        //- Recurrence relation for polynomials orthogonal to kernel function
        virtual void recurrenceRelation
        (
            scalarDiagonalMatrix& a,
            scalarDiagonalMatrix& b,
            scalar primaryAbscissa,
            scalar sigma
        );

        //- Compute maximum value of sigma to ensure realizability
        virtual scalar sigmaMax(univariateMomentSet& moments);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
