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
    Foam::betaEQMOM

Description
    Extended quadrature method of moments with beta kernel density functions.

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
    betaEQMOM.C

\*---------------------------------------------------------------------------*/

#ifndef betaEQMOM_H
#define betaEQMOM_H

#include "extendedMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class betaEQMOM Declaration
\*---------------------------------------------------------------------------*/

class betaEQMOM
:
    public extendedMomentInversion
{
public:

    //- Runtime type information
    TypeName("beta");


    // Constructors

        //- Construct from a dictionary
        betaEQMOM
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    //- Destructor
    virtual ~betaEQMOM();


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
        virtual scalar m2N
        (
            scalar sigma,
            const univariateMomentSet& momentsStar
        );

        //- Recurrence relation for polynomials orthogonal to kernel function
        virtual void recurrenceRelation
        (
            scalarList& a,
            scalarList& b,
            scalar primaryAbscissa,
            scalar sigma
        );

        //- Compute maximum value of sigma to ensure realizability
        virtual scalar sigmaMax(univariateMomentSet& moments);

        //- Return the sum of the kernel density functions at a field of points
        virtual tmp<scalarField> f(const scalarField& x) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
