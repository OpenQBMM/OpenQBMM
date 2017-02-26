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
    Foam::gaussLobattoMomentInversion

Description
    Performs Gauss-Lobatto quadrature on a moment set with an even number of
    moments.

SourceFiles
    gaussLobattoMomentInversion.C
    gaussLobattoMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef gaussLobattoMomentInversion_H
#define gaussLobattoMomentInversion_H

#include "scalarMatrices.H"
#include "univariateMomentInversion.H"
#include "univariateMomentSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class gaussLobattoMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class gaussLobattoMomentInversion
:
    public univariateMomentInversion
{
    // Private data

        //- Bool to indicate the quadrature fell back to Gaussian if a
        //  subset of realizable moments with an even number of elements is
        //  found
        bool forceRadau_;

        //- Minimum known abscissa for Guauss-Lobatto quadrature
        scalar minKnownAbscissa_;

        //- Second known abscissa for Guauss-Lobatto quadrature
        scalar maxKnownAbscissa_;


    // Primate member functions

        //- Update recurrence relation for Gauss-Lobatto quadrature
        void updateRecurrenceRelation();


protected:

    // Protected member functions

        //- Calculate the number of quadrature nodes
        void calcNQuadratureNodes();


public:

    // Constructors

        //- Construct from univariateMomentSet
        gaussLobattoMomentInversion
        (
            univariateMomentSet& moments,
            scalar minKnownAbscissa,
            scalar maxKnownAbscissa
        );

    //- Destructor
    virtual ~gaussLobattoMomentInversion();


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
