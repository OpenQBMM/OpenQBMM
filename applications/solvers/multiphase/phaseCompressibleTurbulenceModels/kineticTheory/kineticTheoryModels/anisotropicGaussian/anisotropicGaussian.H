/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
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
    Foam::anisotropicGaussian

Description


SourceFiles
    anisotropicGaussian.C

\*---------------------------------------------------------------------------*/

#ifndef anisotropicGaussian_H
#define anisotropicGaussian_H

#include "fvMatrix.H"
#include "phaseModel.H"
#include "dragModel.H"
#include "frictionalStressModel.H"
#include "radialModel.H"
#include "AGmomentTransportModel.H"
#include "kineticTheoryModel.H"
#include "fluxSplittingFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{

/*---------------------------------------------------------------------------*\
                           Class anisotropicGaussian Declaration
\*---------------------------------------------------------------------------*/
class anisotropicGaussian
:
    public kineticTheoryModel
{
    // Private data

        dimensionedScalar alphaTheta_;

        dimensionedScalar alphaSigma_;

        //- Modified coefficient of restitution
        dimensionedScalar eta_;

        //- function to calculate flux splitting function
        autoPtr<fluxSplittingFunction> h2Function_;

        //- Switching function
	    volScalarField  h2Fn_;

        //- Particle frictional pressure
        volScalarField PsFric_;

        //- Second order velocity moment tensor
        volSymmTensorField Sigma_;

        //- Moment transport model
		AGmomentTransportModel AGtransport_;


    // Private Member Functions

        //- Update flux weighting function
        void updateh2Fn();


public:

    //- Runtime type information
    TypeName("anisotropicGaussian");

    // Constructors

        //- Construct from components
        anisotropicGaussian
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    ~anisotropicGaussian();

    // Return functions

        //- Return granular stress tensor
        virtual tmp<volSymmTensorField> Sigma() const;


    // Member Functions

        //- Update kinetic theory quantities
        virtual void solve
        (
            const volScalarField& beta,
            const volScalarField& alpha,
            const volTensorField& gradU,
            const volSymmTensorField D
        );

        //- Trasport moments;
        virtual void transportMoments();

        //- Max time-step
        virtual scalar maxUxDx() const;

        //- Return hydrodynamic scale
        virtual tmp<volScalarField> h2() const;

        //- Return interpolated hydrodynamic scale
        virtual tmp<surfaceScalarField> h2f() const;

        //- Return dilute alpha flux
        virtual tmp<volScalarField> ddtAlphaDilute() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
