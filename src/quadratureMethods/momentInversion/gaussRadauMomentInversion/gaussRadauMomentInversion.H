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
    Foam::gaussRadauMomentInversion

Description
    Performs Gauss-Radau quadrature on a moment set with an odd number of
    moments, enforcing that one quadrature abscissa is equal to the specified
    value.

    The quadrature formula degenerates to standard Gaussian quadrature if the
    number of realizable moments is even, in order to preserve the largest
    possible number of moments.

SourceFiles
    gaussRadauMomentInversion.C
    gaussRadauMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef gaussRadauMomentInversion_H
#define gaussRadauMomentInversion_H

#include "scalarMatrices.H"
#include "univariateMomentInversion.H"
#include "univariateMomentSet.H"

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

        //- Bool to indicate the quadrature fell back to Gaussian if a
        //  subset of realizable moments with an even number of elements is
        //  found
        bool forceGauss_;

        //- Known abscissa for Gauss-Radau quadrature
        scalar knownAbscissa_;


    // Primate member functions

        //- Update recurrence relation for Gauss-Radau quadrature
        void updateRecurrenceRelation();


protected:

    // Protected member functions

        //- Calculate the number of quadrature nodes
        void calcNQuadratureNodes();


public:

    // Constructors

        //- Construct from univariateMomentSet
        gaussRadauMomentInversion
        (
            univariateMomentSet& moments,
            scalar knownAbscissa = 0.0
        );

    //- Destructor
    virtual ~gaussRadauMomentInversion();


    // Member functions

        //- Invert moments
        virtual void invert();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
