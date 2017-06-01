/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
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
    Foam::bubblePressureModel

Description

SourceFiles
    bubblePressureModel.C
    newWallLubricationModel.C

\*---------------------------------------------------------------------------*/

#ifndef bubblePressureModel_H
#define bubblePressureModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wallDependentModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                      Class bubblePressureModel Declaration
\*---------------------------------------------------------------------------*/

class bubblePressureModel
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

        //- Coefficient of bubble pressure
        dimensionedScalar Cbp_;


public:

    //- Runtime type information
    TypeName("bubblePressureModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            bubblePressureModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Force dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from components
        bubblePressureModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~bubblePressureModel();


    // Selectors

        static autoPtr<bubblePressureModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Return bubble pressure
        virtual tmp<volScalarField> pb(const label, const label) const = 0;

        //- Return effective viscosity dues to dispersed phase
        virtual tmp<volScalarField> nuEff(const label, const label) const;

        //- Return continuous phase viscous stress
        virtual tmp<volVectorField> divDevRhoReff
        (
            const label,
            const label
        ) const;

        //- Return phase-intensive bubble pressure force
        virtual tmp<volVectorField> Fi(const label, const label) const;

        //- Return bubble pressure force
        virtual tmp<volVectorField> F(const label, const label) const;

        //- Return face bubble pressure force
        virtual tmp<surfaceScalarField> Ff(const label, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
