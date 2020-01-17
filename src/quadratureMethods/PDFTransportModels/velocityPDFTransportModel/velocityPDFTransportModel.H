/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2018 by Alberto Passalacqua
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
    Foam::velocityPDFTransportModel

Description
    Solve a velocity PDF transport equation with the hyperbolic conditional
    quadrature method of moments

SourceFiles
    velocityPDFTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef velocityPDFTransportModel_H
#define velocityPDFTransportModel_H

#include "PDFTransportModel.H"
#include "quadratureApproximations.H"
#include "velocityMomentAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{

/*---------------------------------------------------------------------------*\
                Class velocityPDFTransportModel Declaration
\*---------------------------------------------------------------------------*/

class velocityPDFTransportModel
:
    public PDFTransportModel
{
protected:

    // Protected data

        //- velocity quadrature approximation used to solve PDF transport
        velocityQuadratureApproximation quadrature_;

        //- Moment advection
        autoPtr<velocityMomentAdvection> momentAdvection_;

    // Protected member functions

        //- Update implicit collison source terms
        virtual void updateImplicitMomentSource() = 0;

        //- Calculate implicit source terms
        virtual tmp<fvScalarMatrix> implicitMomentSource
        (
            const volVelocityMoment& moment
        ) = 0;

        //- Calculate explicit source terms with realizable ODE solver
        virtual void explicitMomentSource() = 0;

        //- Are moments updated due to sources
        virtual bool solveMomentSources() const = 0;

        //- Is the realizable ode solver used to solve sources
        virtual bool solveMomentOde() const = 0;


public:

    // Constructors

        //- Construct from components
        velocityPDFTransportModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const word& support
        );


    //- Destructor
    virtual ~velocityPDFTransportModel();


    // Member Functions

        //- Access

            //- Returns a const reference to the quadrature approximation
            virtual inline const velocityQuadratureApproximation& quadrature() const;

            //- Return the maximum Courant number ensuring moment realizability
            virtual inline scalar realizableCo() const;

            //- Return the maximum Courant number based on velocities
            virtual inline scalar CoNum() const;


        //-Edit

            //- Solve velocity PDF transport equation
            virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PDFTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "velocityPDFTransportModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
