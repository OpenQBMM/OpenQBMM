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
    Foam::hyperbolicMomentInversion

Description
    Performs analytical HyQMOM inversion on a moment set

SourceFiles
    hyperbolicMomentInversion.C
    hyperbolicMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef hyperbolicMomentInversion_H
#define hyperbolicMomentInversion_H

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
    //- Central moments
    scalarList centralMoments_;

    //- Calculate scaled and central moments
    void calculateCentralMoments(const scalarList& moments);

    //- Scale central abscissae to normal
    void unscaleAbscissae(const scalarList& moments);


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

        //- Calculates the number of quadrature nodes
        virtual void calcNQuadratureNodes
        (
            univariateMomentSet& moments
        );

        //- Invert the set of moments to compute weights and abscissae
        virtual void invert
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
