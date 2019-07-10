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
    Foam::fv::phaseIncompressibleMeanVelocityForce

Description
    Creates a pressure gradient source

    Note: Currently only handles kinematic pressure

    \heading Source usage
    Example usage:
    \verbatim
    \verbatim
    phaseIncompressibleMeanVelocityForceCoeffs
    {
        selectionMode    all;                    // Apply force to all cells
        fieldNames       (U.gas);                //
        alphaName         alpha.gas;
        Ubar             (10.0 0 0);             // Desired mean velocity
        relaxation       0.2;                    // Optional relaxation factor
    }
    \endverbatim

SourceFiles
    phaseIncompressibleMeanVelocityForce.C

\*---------------------------------------------------------------------------*/

#ifndef phaseIncompressibleMeanVelocityForce_H
#define phaseIncompressibleMeanVelocityForce_H

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
               Class phaseIncompressibleMeanVelocityForce Declaration
\*---------------------------------------------------------------------------*/

class phaseIncompressibleMeanVelocityForce
:
    public cellSetOption
{
    // Private data
        // alpha name
        const volScalarField& alpha_;

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
    TypeName("phaseIncompressibleMeanVelocityForce");


    // Constructors

        //- Construct from explicit source name and mesh
        phaseIncompressibleMeanVelocityForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        phaseIncompressibleMeanVelocityForce
        (
            const phaseIncompressibleMeanVelocityForce&
        ) = delete;


    // Member Functions

        // Evaluate

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
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
        void operator=(const phaseIncompressibleMeanVelocityForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
