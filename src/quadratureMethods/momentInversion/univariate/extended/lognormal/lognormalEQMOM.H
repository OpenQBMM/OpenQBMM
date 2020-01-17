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
    Foam::lognormalEQMOM

Description
    Extended quadrature method of moments with log-normal kernel density
    functions.

    References
    \verbatim
        "An extended quadrature-based moment method with log-normal kernel
        density functions"
        E Madadi-Kandjani, A Passalacqua
        Chemical Engineering Science
        Volume 131, Pages 323-339, 2015
    \endverbatim

    \verbatim
        "The theory of canonical moments with applications in Statistics,
        Probability and Analysis"
        H Dette, W J Studden
        Wiley & Sons, 1997
    \endverbatim

    \verbatim
        "A general approximation method for solving integrals containing a
        lognormal weighting function"
        M Wilck
        Journal of Aerosol Science
        Volume 32, Pages 1111-1116, 2001
    \endverbatim

SourceFiles
    lognormalEQMOM.C

\*---------------------------------------------------------------------------*/

#ifndef lognormalEQMOM_H
#define lognormalEQMOM_H

#include "extendedMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class lognormalEQMOM Declaration
\*---------------------------------------------------------------------------*/

class lognormalEQMOM
:
    public extendedMomentInversion
{
public:

    //- Runtime type information
    TypeName("lognormal");


    // Constructors

        //- Construct from a dictionary and label
        lognormalEQMOM
        (
            const dictionary& dict,
            const label nMoments,
            const label nSecondaryNodes
        );


    //- Destructor
    virtual ~lognormalEQMOM();


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
