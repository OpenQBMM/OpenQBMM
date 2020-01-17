/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::univariateMomentAdvection

Description
    Abstract class for the advection of a univariate moment set.

SourceFiles
    univariateMomentAdvection.C
    univariateMomentAdvectionI.H

\*---------------------------------------------------------------------------*/

#ifndef univariateMomentAdvection_H
#define univariateMomentAdvection_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "mappedPtrList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "quadratureNodes.H"
#include "momentFieldSets.H"
#include "basicFieldMomentInversion.H"
#include "quadratureApproximations.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class univariateMomentAdvection Declaration
\*---------------------------------------------------------------------------*/

class univariateMomentAdvection
{
protected:

    // Protected data

        //- Name of the distribution associated to the quadrature approximation.
        const word name_;

        //- Reference to the moments to advect
        const volScalarMomentFieldSet& moments_;

        //- Number of moments
        const label nMoments_;

        //- Advection terms for moment transport equations
        mappedPtrList<volScalarField> divMoments_;

        //- Field used to find upwind values of the owner
        surfaceScalarField own_;

        //- Field used to find upwind values of the neighbour
        surfaceScalarField nei_;

        //- Conservative advection flux
        const surfaceScalarField& phi_;

        //- Support of the distribution function
        word support_;

        //- Dimensions of the underlying distribution
        label nDimensions_;

        //- Construct field name from name, order
        word fieldName
        (
            const word& name,
            const labelList& order
        ) const;


public:

    //- Runtime type information
    TypeName("univariateMomentAdvection");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMomentAdvection,
            dictionary,
            (
                const dictionary& dict,
                const scalarQuadratureApproximation& quadrature,
                const surfaceScalarField& phi,
                const word& support
            ),
            (dict, quadrature, phi, support)
        );


    // Constructors

        //- Construct from univariateMomentSet
        univariateMomentAdvection
        (
            const dictionary& dict,
            const scalarQuadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );

        //- Disallow default bitwise copy construct
        univariateMomentAdvection(const univariateMomentAdvection&) = delete;


    //- Destructor
    virtual ~univariateMomentAdvection();


    // Selectors

        static autoPtr<univariateMomentAdvection> New
        (
            const dictionary& dict,
            const scalarQuadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    // Member Functions

        //- Return the advection term for the moment transport equation
        inline const mappedPtrList<volScalarField>& divMoments() const;

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const = 0;

        //- Update moment advection
        virtual void update() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const univariateMomentAdvection&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentAdvectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
