/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
    Foam::pdPhaseModel

Description
    Class for a polydisperse phase model. Each bubble size hase a unique mean
    velocity. Size and velocity moments are stored in quadrature.

SourceFiles
    pdPhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef pdPhaseModel_H
#define pdPhaseModel_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "phaseModel.H"
#include "monoKineticQuadratureApproximation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class pdPhaseModel
:
    public phaseModel
{
    // Private data

        //- QuadratureApproximation
        monoKineticQuadratureApproximation quadrature_;

        //- Number of nodes
        const label nNodes_;

        //- Number of moments
        const label nMoments_;

        //- List of alphas
        PtrList<volScalarField> alphas_;

        //- List of velocities (Refrence to the velocities stored in the
        //  quadrature
        PtrList<volVectorField>& Us_;

        //- List of deviation velocities
        PtrList<volVectorField> Vs_;

        //- List of diameters
        PtrList<volScalarField> ds_;

        //- Maximum diameter
        dimensionedScalar maxD_;

        //- Minimum diameter
        dimensionedScalar minD_;


public:

    // Constructors
        pdPhaseModel
        (
            const fvMesh& mesh,
            const dictionary& phaseProperties,
            const word& phaseName
        );


    //- Destructor
    virtual ~pdPhaseModel();


    // Member Functions

        //- Return the number of nodes
        virtual label nNodes() const
        {
            return nNodes_;
        }

        //- Return the diameter for nodei
        virtual const volScalarField& ds(const label nodei) const
        {
            return ds_[nodei];
        }

        //- Return alpha field for nodei
        virtual const volScalarField& alphas(const label nodei) const
        {
            return alphas_[nodei];
        }

        //- Return non-constant access to alpha field for nodei
        virtual volScalarField& alphas(const label nodei)
        {
            return alphas_[nodei];
        }

        //- Return the velocity for nodei
        virtual const volVectorField& Us(const label nodei) const
        {
            return Us_[nodei];
        }

        //- Return non-const access to the velocity for nodei
        virtual volVectorField& Us(const label nodei)
        {
            return Us_[nodei];
        }

        //- Return deviation velocity for nodei
        virtual tmp<volVectorField> Vs(const label nodei)
        {
            return Vs_[nodei];
        }

        //- Correct the phase properties
        virtual void correct();

        //- Relative transport of moments
        virtual void relativeTransport();

        //- Average transport of moments
        virtual void averageTransport(const PtrList<fvVectorMatrix>& AEqns);

        //- Update moments
        virtual void updateMoments();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
