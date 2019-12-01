/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::fv::phaseCompressibleMeanVelocityForce

Description
    Calculates and applies the force necessary to maintain the specified mean
    velocity.

    Note: Currently only handles kinematic pressure (incompressible solvers).

    \heading Source usage
    Example usage:
    \verbatim
    phaseCompressibleMeanVelocityForceCoeffs
    {
        selectionMode    all;                    // Apply force to all cells
        fieldNames       (U.gas);                //
        alphaName         alpha.gas;
        rhoName           rho.gas;
        Ubar             (10.0 0 0);             // Desired mean velocity
        relaxation       0.2;                    // Optional relaxation factor
    }
    \endverbatim

SourceFiles
    phaseCompressibleMeanVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef phaseCompressibleMeanVelocityForce_H
#define phaseCompressibleMeanVelocityForce_H

#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class phaseCompressibleMeanVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class phaseCompressibleMeanVelocityForce
:
    public cellSetOption
{
    // Private data

        const volScalarField& alpha_;

        const volScalarField& rho_;

        //- Average velocity
        vector Ubar_;

        //- target velocity magnitude
        scalar magUbar_;

        //- Flow direction
        vector flowDir_;

        //- Relaxation factor
        scalar relaxation_;

        //- Pressure gradient before correction
        scalar gradP0_;

        //- Change in pressure gradient
        scalar dGradP_;

        //- Matrix 1/A coefficients field pointer
        autoPtr<volScalarField> rAPtr_;

    // Private Member Functions

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);


public:

    //- Runtime type information
    TypeName("phaseCompressibleMeanVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        phaseCompressibleMeanVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        phaseCompressibleMeanVelocityForce
        (
            const phaseCompressibleMeanVelocityForce&
        )= delete;


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldI
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldI
            );


        // I-O

            //- Read source dictionary
            virtual bool read(const dictionary& dict)
            {
                return false;
            }

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseCompressibleMeanVelocityForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
