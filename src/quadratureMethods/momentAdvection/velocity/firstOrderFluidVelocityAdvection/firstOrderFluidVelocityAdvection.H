/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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
    Foam::velocityAdvection::firstOrderFluidVelocity

Description
    Transport of the moments of a distribution whose nodes all move with one
    velocity, the velocity of the fluid they are carried by, as the mixture
    fraction and the compositions conditioned on it are.

    Every moment then obeys the same conservation law, d m / d t + div(U m)
    = 0, and there is no flux to split by node: the moments of a face are
    the moments of the cell upwind of it. Taking them from the upwind cell
    is what the kinetic schemes do when every node has the same velocity,
    and it leaves the moment set of a cell a convex combination of the sets
    it was built from, which is realizable, up to a Courant number of one.
    The conditional means are convex combinations of those of the upwind
    cells for the same reason, so they stay within their support.

    The flux is the one of the fluid, named by phi in the dictionary of the
    advection, phi unless another name is given.

    Where that flux is a mass flux, which it is in a solver of variable
    density, the distribution is mass-weighted: its weights are fractions
    of the mass of the fluid, and what is conserved is rho m,

        d (rho m) / d t + div(rho U m) = 0.

    The scheme then gives div(rho U m), upwind as before, and a transport
    model balances it with the rate of change of rho m. With the density of
    a cell obeying the same discrete continuity equation with the same
    flux, which is what a solver obtains by solving it before the moments,
    the new moments of a cell are again a convex combination of those of
    the cell and of the cells upwind of it, up to a Courant number of one
    taken with the mass flux over the mass of the cell, and the moment of
    order zero of a distribution that integrates to one stays one. The
    density is named by rho, rho unless another name is given.

    firstOrderKinetic is the scheme for a distribution whose nodes have
    velocities of their own, which it splits the flux of a face by; this one
    is the same reconstruction where they share the velocity of the fluid,
    and is what the univariate firstOrderKinetic does for a quadrature of
    one dimension.

    This is a first-order scheme.

SourceFiles
    firstOrderFluidVelocityAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef firstOrderFluidVelocityAdvection_H
#define firstOrderFluidVelocityAdvection_H

#include "velocityMomentAdvection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{

/*---------------------------------------------------------------------------*\
                       Class firstOrderFluidVelocity Declaration
\*---------------------------------------------------------------------------*/

class firstOrderFluidVelocity
:
    public velocityMomentAdvection
{
    // Private data

        //- Name of the flux of the fluid that carries the distribution
        const word phiName_;

        //- Name of the density of the fluid, read where the flux is a
        //  mass flux
        const word rhoName_;

        //- The mesh, to find the flux on
        const fvMesh& mesh_;


    // Private Member Functions

        //- The flux of the fluid
        const surfaceScalarField& phi() const;


public:

    //- Runtime type information
    TypeName("firstOrderFluidVelocity");


    // Constructors

        //- Construct from a dictionary and a quadrature approximation
        firstOrderFluidVelocity
        (
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            const List<supportType>& supports
        );


    //- Destructor
    virtual ~firstOrderFluidVelocity();


    // Member Functions

        //- Whether the flux of the fluid is a mass flux
        virtual bool densityWeighted() const;

        //- The largest Courant number the moments stay realizable at: the
        //  moments of a cell are a convex combination of those it is
        //  updated from up to one
        virtual scalar realizableCo() const;

        //- The Courant number of the flux of the fluid
        virtual scalar CoNum() const;

        //- Update the advection terms of the moments
        virtual void update();

        //- Update them with the flux given instead of the one of the fluid
        virtual void update
        (
            const surfaceScalarField& phi,
            const bool wallCollisions
        );

        //- Not used by a distribution whose nodes share their velocity
        virtual void update
        (
            const mappedPtrList<volVectorField>& Us,
            const bool wallCollisions
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace velocityAdvection
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
