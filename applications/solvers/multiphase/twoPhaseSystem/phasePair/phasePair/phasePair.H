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
    Foam::phasePair

Description

SourceFiles
    phasePair.C

\*---------------------------------------------------------------------------*/

#ifndef phasePair_H
#define phasePair_H

#include "phaseModel.H"
#include "phasePairKey.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class phasePair Declaration
\*---------------------------------------------------------------------------*/

class phasePair
:
    public phasePairKey
{
public:

    // Hash table types

        //- Dictionary hash table
        typedef HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

        //- Scalar hash table
        typedef HashTable<scalar, phasePairKey, phasePairKey::hash>
            scalarTable;


private:

    // Private data

        //- Phase 1
        const phaseModel& phase1_;

        //- Phase 2
        const phaseModel& phase2_;

        //- Gravitational acceleration
        const dimensionedVector& g_;

        //- Surface tension coefficient
        const dimensionedScalar sigma_;


    // Private member functions

        // Etvos number for given diameter
        tmp<volScalarField> EoH(const volScalarField& d) const;

public:

    // Constructors

        //- Construct from two phases, gravity and surface tension table
        phasePair
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const dimensionedVector& g,
            const scalarTable& sigmaTable,
            const bool ordered = false
        );


    //- Destructor
    virtual ~phasePair();


    // Member Functions

        //- Dispersed phase
        virtual const phaseModel& dispersed() const;

        //- Continuous phase
        virtual const phaseModel& continuous() const;

        //- Pair name
        virtual word name() const;

        //- Average density
        tmp<volScalarField> rho() const;

        //- Relative velocity magnitude
        tmp<volScalarField> magUr(const label nodei, const label nodej) const;

        //- Relative velocity magnitude
        tmp<volScalarField> magUr() const;

        //- Relative velocity
        tmp<volVectorField> Ur(const label nodei, const label nodej) const;

        //- Mean relative velocity
        tmp<volVectorField> Ur() const;

        //- Reynolds number
        tmp<volScalarField> Re(const label nodei, const label nodej) const;

        //- Mean Reynolds number
        tmp<volScalarField> Re() const;

        //- Weber number
        tmp<volScalarField> We(const label nodei, const label nodej) const;

        //- Mean Weber number
        tmp<volScalarField> We() const;

        //- Eotvos number
        tmp<volScalarField> Eo(const label nodei, const label nodej) const;

        //- Mean Eotvos number
        tmp<volScalarField> Eo() const;

        //- Eotvos number based on hydraulic diameter type 1
        tmp<volScalarField> EoH1(const label nodei, const label nodej) const;

        //- Mean Eotvos number based on hydraulic diameter type 1
        tmp<volScalarField> EoH1() const;

        //- Eotvos number based on hydraulic diameter type 2
        tmp<volScalarField> EoH2(const label nodei, const label nodej) const;

        //- Mean Eotvos number based on hydraulic diameter type 2
        tmp<volScalarField> EoH2() const;

        //- Morton Number
        tmp<volScalarField> Mo() const;

        //- Takahashi Number
        tmp<volScalarField> Ta(const label nodei, const label nodej) const;

        //- Mean Takahashi Number
        tmp<volScalarField> Ta() const;

        //- Aspect ratio
        virtual tmp<volScalarField> E
        (
            const label nodei,
            const label nodej
        ) const;

        //- Mean Aspect ratio
        virtual tmp<volScalarField> E() const;

        // Access

            // Phase 1
            inline const phaseModel& phase1() const;

            // Phase 2
            inline const phaseModel& phase2() const;

            // Gravitational acceleration
            inline const dimensionedVector& g() const;

            // Surface tension coefficient
            inline const dimensionedScalar& sigma() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phasePairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
