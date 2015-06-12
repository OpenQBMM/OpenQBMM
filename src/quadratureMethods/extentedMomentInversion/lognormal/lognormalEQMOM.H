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
            const univariateMomentSet& moments,
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
        univariateMomentSet& moments
    );

    //- Compute the starred moments
    virtual void momentsToMomentsStar(scalar sigma);

    //- Recurrence relation for polynomials orthogonal to kernel function
    virtual void recurrenceRelation
    (
        scalarDiagonalMatrix& a, 
        scalarDiagonalMatrix& b,
        scalar primaryAbscissa
    );
    
    //- Compute maximum value of sigma to ensure realizability
    virtual scalar sigmaMax();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
