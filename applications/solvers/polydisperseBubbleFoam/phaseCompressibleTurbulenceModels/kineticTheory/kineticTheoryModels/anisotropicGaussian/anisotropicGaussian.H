/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::kineticTheoryModel

Description


SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "fvMatrix.H"
#include "phaseModel.H"
#include "dragModel.H"
#include "frictionalStressModel.H"
#include "radialModel.H"
#include "anisotropicGaussianModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
    :
    public IOdictionary
{
        // Private data

        // Input Fields
        const fvMesh& mesh_;

        phaseModel& particles_;

        const phaseModel& gas_;
        
        // Kinetic Theory Model coefficients

        //- Coefficient of restitution
        dimensionedScalar e_;
        
        dimensionedScalar eta_;

        //- Maximum packing phase-fraction
        dimensionedScalar alphaMax_;

        //- Min value for which the frictional stresses are zero
        dimensionedScalar alphaMinFriction_;
        
        //- Residual phase fraction
        dimensionedScalar residualAlpha_;
        
        // Sub-models
        //- Run-time selected frictional stress model
        autoPtr<kineticTheoryModels::frictionalStressModel>
        frictionalStressModel_;
       
        //- Run-time selected radial distribution model
        autoPtr<kineticTheoryModels::radialModel> radialModel_;
        
        
        // Kinetic Theory Model Fields
	
        const volScalarField& alphap_;
        
        //- The granular energy/temperature
        volScalarField Theta_;

		//- The particle viscosity
        volScalarField nut_;   
        
        //- The granular bulk viscosity
        volScalarField lambda_;

        //- The granular radial distribution
        volScalarField gs0_;

        //- The granular "thermal" conductivity
        volScalarField kappa_;
        
        // Frictional pressure
        volScalarField ppfr_;
		
	    volScalarField&  h2Fn_;
		        
        volSymmTensorField Sigma_;
		
		anisotropicGaussianModel diluteAGmodel_;
		
		const word  h2FnMethod_;
		
		const label h2FnParaPow_;
				
        // Private Member Functions

        //- Disallow default bitwise copy construct
        kineticTheoryModel(const kineticTheoryModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&);


    public:


        //- Construct from components
        kineticTheoryModel
        (
            const fvMesh& mesh,
            phaseModel& particles,
            const phaseModel& gas,
            volScalarField&  h2Fn
        );


        //- Destructor
        ~kineticTheoryModel();

        // Member Functions

        //- Re-read model coefficients if they have changed
        bool read();
        
		tmp<volScalarField> pPrime() const;
		
        //- Return the effective stress tensor
        tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        void updateh2Fn() ;
        
        void updateViscosity(const volScalarField& K);
        
        //- Solve the kinetic theory equations and correct the viscosity
        void solveDense
        (
			const volScalarField& K,
			const volScalarField& particleContinuityErr
		);
        
        void solveDilute
        (
			const surfaceScalarField& h2f
        );
		
        volScalarField& Theta()
        {
            return Theta_;
        }
		
		volSymmTensorField& Sigma()
        {
            return Sigma_;
        }
        
        const Switch& AGmodel()
        {
            return diluteAGmodel_.AGmodel();
        };

        const scalar& alphaMax() const
        {
            return alphaMax_.value();
        }
        
        scalar maxUxDx() 
        {
			return diluteAGmodel_.maxUxDx();
		}


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
