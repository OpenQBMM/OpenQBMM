/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added additional return functions so that class can
                            be extended to polydisperse
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
    Foam::phaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{
protected:
    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Name of phase
        word name_;

        //- Phase dictionary
        dictionary phaseDict_;

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        dimensionedScalar residualAlpha_;

        //- Density
        dimensionedScalar rho_;

        //- Viscosity
        dimensionedScalar nu_;

        //- Mean velocity
        volVectorField U_;

        //- Volumetric flux of the phase
        surfaceScalarField alphaPhi_;

        //- Mass flux of the phase
        surfaceScalarField alphaRhoPhi_;

        //- Volumetric flux of the phase
        autoPtr<surfaceScalarField> phiPtr_;

        //- Diameter (alpha weighted diameter for polydisperse)
        volScalarField d_;


public:

    // Constructors
        phaseModel
        (
            const fvMesh& mesh,
            const dictionary& phaseProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~phaseModel();


    // Member Functions

        //- Return the name of this phase
        const word& name() const
        {
            return name_;
        }

        //- Return the number of nodes
        virtual label nNodes() const
        {
            return 1;
        }

        //- Return the residual phase-fraction for given phase
        //  Used to stabilize the phase momentum as the phase-fraction -> 0
        const dimensionedScalar& residualAlpha() const
        {
            return residualAlpha_;
        }

        //- Return the mean diameter
        virtual const volScalarField& d() const
        {
            return d_;
        }

        //- Return diameter of nodei
        virtual const volScalarField& ds(const label nodei) const
        {
            return d_;
        }

        //- Return the density
        tmp<volScalarField> rho() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "rho",
                            name_
                        ),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_,
                    rho_
                )
            );
        }

        //- Return the viscosity
        tmp<volScalarField> nu() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "nu",
                            name_
                        ),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_,
                    nu_
                )
            );
        }

        //- Return the viscosity
        tmp<volScalarField> mu() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "mu",
                            name_
                        ),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_,
                    nu_*rho_
                )
            );
        }

        //- Constant access to alpha field for nodei
        virtual const volScalarField& alphas(const label nodei) const
        {
            return *this;
        }

        //- Non-constant access to alpha field for nodei
        virtual volScalarField& alphas(const label nodei)
        {
            return *this;
        }

        //- Constant access to the mean velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Non-const access to the mean velocity
        volVectorField& U()
        {
            return U_;
        }

        //- Constant access to the velocity for nodei
        virtual const volVectorField& Us(const label nodei) const
        {
            return U_;
        }

        //- Non-const access to the velocity for nodei
        virtual volVectorField& Us(const label nodei)
        {
            return U_;
        }

        //- Deviation of velocity nodei from the mean
        virtual tmp<volVectorField> Vs(const label nodei)
        {
            return tmp<volVectorField>
            (
                new volVectorField
                (
                    IOobject
                    (
                        IOobject::groupName
                        (
                            "V",
                            name_
                        ),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_,
                    dimensionedVector("Zero",dimVelocity,Zero)
                )
            );
        }

        //- Constant access to the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phiPtr_();
        }

        //- Non-const access to the volumetric flux
        surfaceScalarField& phi()
        {
            return phiPtr_();
        }

        //- Return the volumetric flux of the phase
        const surfaceScalarField& alphaPhi() const
        {
            return alphaPhi_;
        }

        //- Return non-const access to the volumetric flux of the phase
        surfaceScalarField& alphaPhi()
        {
            return alphaPhi_;
        }

        //- Return the mass flux of the phase
        const surfaceScalarField& alphaRhoPhi() const
        {
            return alphaRhoPhi_;
        }

        //- Return non-const access to the mass flux of the phase
        surfaceScalarField& alphaRhoPhi()
        {
            return alphaRhoPhi_;
        }

        //- Ensure that the flux at inflow BCs is preserved
        void correctInflowFlux(surfaceScalarField& alphaPhi) const;

        //- Correct the phase properties
        virtual void correct();

        //- Relative transport of moments, not used in standard phase model
        virtual void relativeTransport()
        {
            return;
        }

        //- Average transport of moments, not used in standard phase model
        virtual void averageTransport(PtrList<fvVectorMatrix>& AEqns)
        {
            return;
        }

        virtual void updateMoments()
        {
            return;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
