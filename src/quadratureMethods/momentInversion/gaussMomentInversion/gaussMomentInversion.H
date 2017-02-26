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
    Foam::gaussMomentInversion

Description
    Performs Gaussian quadrature on a moment set with an even numer of moments.

SourceFiles
    gaussMomentInversion.C
    gaussMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef gaussMomentInversion_H
#define gaussMomentInversion_H

#include "scalarMatrices.H"
#include "univariateMomentInversion.H"
#include "univariateMomentSet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class gaussMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class gaussMomentInversion
:
    public univariateMomentInversion
{
protected:

    // Protected member functions

        //- Calculate the number of quadrature nodes
        void calcNQuadratureNodes();


public:

    // Constructors

        //- Construct from univariateMomentSet
        gaussMomentInversion
        (
            univariateMomentSet& moments
        );

    //- Destructor
    virtual ~gaussMomentInversion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
