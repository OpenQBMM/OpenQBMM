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

        dimensionedScalar eta_;

        // Frictional pressure
        volScalarField ppfr_;

	    volScalarField  h2Fn_;

        volSymmTensorField Sigma_;

		AGmomentTransportModel AGtransport_;

		const word  h2FnMethod_;

		const scalar h2FnParaPow_;

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

    // Member Functions

        //- Update viscosities
        virtual void updateViscosities();

        //- Update kinetic theory quantities
        virtual void correct();

        //- Trasport moments;
        virtual void transportMoments();

        //- Max time-step
        virtual scalar maxUxDx() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticTheoryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
