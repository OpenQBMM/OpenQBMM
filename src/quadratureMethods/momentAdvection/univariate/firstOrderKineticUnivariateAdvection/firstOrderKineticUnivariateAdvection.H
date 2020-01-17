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
    Foam::firstOrderKinetic

Description
    First-order univariate moment advection with kinetic scheme.

SourceFiles
    firstOrderKineticUnivariateAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef firstOrderKineticUnivariateAdvection_H
#define firstOrderKineticUnivariateAdvection_H

#include "univariateMomentAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace univariateAdvection
{

/*---------------------------------------------------------------------------*\
                    Class firstOrderKinetic Declaration
\*---------------------------------------------------------------------------*/

class firstOrderKinetic
:
    public univariateMomentAdvection
{
    // Private data

        //- List of quadrature nodes used for moment advection
        autoPtr<mappedPtrList<volScalarNode>> nodes_;

        //- List of interpolated nodes (neighbour)
        autoPtr<mappedPtrList<surfaceScalarNode>> nodesNei_;

        //- List of interpolated nodes (owner)
        autoPtr<mappedPtrList<surfaceScalarNode>> nodesOwn_;

        //- List of interpolated moments (neighbour)
        surfaceScalarMomentFieldSet momentsNei_;

        //- List of interpolated moments (owner)
        surfaceScalarMomentFieldSet momentsOwn_;

        //- Flux moment inversion method
        autoPtr<basicFieldMomentInversion> momentFieldInverter_;

        //- Number of nodes used for advection
        label nNodes_;


protected:

        //- Protected member functions
        void interpolateNodes();


public:

    //- Runtime type information
    TypeName("firstOrderKinetic");


    // Constructors

        //- Construct from univariateMomentSet
        firstOrderKinetic
        (
            const dictionary& dict,
            const scalarQuadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~firstOrderKinetic();


    // Public member functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo() const;

        //- Update moment advection
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace univariateAdvection
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
