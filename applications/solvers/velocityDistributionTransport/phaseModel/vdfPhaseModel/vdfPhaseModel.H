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
    Foam::vdfPhaseModel

Description
    Base class for a moving phase model. Functions are made so that the class
    can be abstracted to a polydisperse phase.

SourceFiles
    vdfPhaseModel.C
    newvdfPhaseModel.C
    vdfPhaseModels.C

\*---------------------------------------------------------------------------*/

#ifndef vdfPhaseModel_H
#define vdfPhaseModel_H

#include "phaseModel.H"
#include "populationBalanceModel.H"
#include "quadratureApproximations.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class vdfPhaseModel Declaration
\*---------------------------------------------------------------------------*/

class vdfPhaseModel
:
    public phaseModel
{
protected:
    // Protected data

        //- Population balance dictionary
        IOdictionary pbeDict_;

        //- Population balance model
        autoPtr<populationBalanceModel> populationBalance_;

        //- Reference to quadrature
        velocityQuadratureApproximation& quadrature_;

        //- PtrList of volume fractions
        PtrList<volScalarField> volumeFractions_;

        //- Are sigma and theta computed
        bool computeVariance_;

        //- Velocity covariance tensor
        tmp<volTensorField> sigma_;

        //- Granular temperature
        tmp<volScalarField> Theta_;

        //- PtrList of diameters
        PtrList<volScalarField> ds_;

        //- Minimum diameter
        dimensionedScalar minD_;


public:

    // Constructors
        vdfPhaseModel
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName
        );

         //- Return clone
        autoPtr<vdfPhaseModel> clone() const;


    //- Destructor
    virtual ~vdfPhaseModel();


    // Member Functions

        //- Return the number of nodes
        virtual label nNodes() const;

        //- Constant access to alpha field for nodei
        virtual tmp<volScalarField> volumeFraction(const label nodei = -1) const;

        //- Return diameter of nodei
        virtual tmp<volScalarField> d(const label nodei = -1) const;

        //- Constant access to the mean velocity
        virtual const volVectorField& U(const label nodei = -1) const;

        //- Non-const access to the mean velocity
        virtual volVectorField& U(const label nodei = -1);

        //- Deviation of velocity nodei from the mean
        virtual tmp<volVectorField> Vs(const label nodei = -1) const;

        //- Solve the underlying population balance model
        virtual void solve();

        //- Correct mean quantities
        virtual void correct();

        //- Relative transport of moments, not used in standard phase model
        virtual void relativeTransport();

        //- Average transport of moments, not used in standard phase model
        virtual void averageTransport();

        //- Solve change in size moments based on breakup and coalesence
        virtual void solveSource();

        //- Return the maximum realizable courant number
        virtual scalar realizableCo() const;

        //- Return the max courant number
        virtual scalar CoNum() const;

        //- Read phase properties dictionary
        virtual bool read()
        {
            return true;
        }

        //- Read phase properties dictionary
        virtual bool read(const bool readOK)
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
