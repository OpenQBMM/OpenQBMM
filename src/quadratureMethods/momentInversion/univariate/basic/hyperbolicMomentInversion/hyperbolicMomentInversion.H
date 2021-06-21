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
    Foam::hyperbolicMomentInversion

Description
    Performs three-node hyperbolic moment inversion given five moments.

SourceFiles
    hyperbolicMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef gaussMomentInversion_H
#define gaussMomentInversion_H

#include "univariateMomentInversion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class hyperbolicMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class hyperbolicMomentInversion
:
    public univariateMomentInversion
{
    // Private data

        //- Minimum value of the variance to compute eta
        scalar etaMin_;

        //- Maximum normalized skewness
        scalar qMax_;

        //- Small negative value for realizability
        scalar smallNegRealizability_;


public:

    //- Runtime type information
    TypeName("hyperbolic");


    // Constructors

        //- Construct from univariateMomentSet
        hyperbolicMomentInversion
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~hyperbolicMomentInversion();


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

        //- Invert the set of moments to compute weights and abscissae
        virtual void invert
        (
            univariateMomentSet& moments,
            const scalar minKnownAbscissa = 0,
            const scalar maxKnownAbscissa = 0
        );

        //- Calculates the number of quadrature nodes
        virtual void calcNQuadratureNodes
        (
            univariateMomentSet& moments
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
