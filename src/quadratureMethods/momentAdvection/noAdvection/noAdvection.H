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
    Foam::noAdvection

Description
    Empty scheme for cases without moment advection.

SourceFiles
    noAdvection.C
    noAdvectionI.H

\*---------------------------------------------------------------------------*/

#ifndef noAdvection_H
#define noAdvection_H

#include "univariateMomentAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class noAdvection Declaration
\*---------------------------------------------------------------------------*/

class noAdvection
:
    public univariateMomentAdvection
{
public:

    //- Runtime type information
    TypeName("noAdvection");


    // Constructors

        //- Construct from univariateMomentSet
        noAdvection
        (
            const dictionary& dict,
            const univariateQuadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~noAdvection();


    // Public member functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo();

        //- Update moment advection
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
