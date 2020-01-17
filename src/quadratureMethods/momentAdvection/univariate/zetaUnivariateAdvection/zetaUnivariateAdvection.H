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
    Foam::zeta

Description
    Second-order univariate moment advection with zeta kinetic scheme.

    References
    \verbatim
        "Realizable second-order finite-volume schemes for the advection of
        moment sets of the particle size distribution"
        F Laurent, T. T. Nguyen
        Journal of Computational Physics
        Volume 337, Pages 309-338, 2017
    \endverbatim

SourceFiles
    zetaUnivariateAdvection.C

\*---------------------------------------------------------------------------*/

#ifndef zetaUnivariateAdvection_H
#define zetaUnivariateAdvection_H

#include "univariateMomentAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace univariateAdvection
{

/*---------------------------------------------------------------------------*\
                            Class zeta Declaration
\*---------------------------------------------------------------------------*/

class zeta
:
    public univariateMomentAdvection
{
    // Private data

        //- Reference to zero-order moment field
        const volScalarField& m0_;

        //- Reconstructed m0 (owner)
        surfaceScalarField m0Own_;

        //- Reconstructed m0 (neighbour)
        surfaceScalarField m0Nei_;

        //- Number of zeta_k values
        label nZetas_;

        //- List of fields of zeta_k (n fields for n + 1 moments)
        PtrList<volScalarField> zetas_;

        //- List of interpolated zeta_k values (neighbour)
        PtrList<surfaceScalarField> zetasNei_;

        //- List of interpolated nodes (owner)
        PtrList<surfaceScalarField> zetasOwn_;

        //- List of interpolated zeta_k values (neighbour)
        PtrList<surfaceScalarField> zetasUpwindNei_;

        //- List of interpolated nodes (owner)
        PtrList<surfaceScalarField> zetasUpwindOwn_;

        //- List of interpolated zeta_k values (neighbour)
        PtrList<surfaceScalarField> zetasCorrNei_;

        //- List of interpolated nodes (owner)
        PtrList<surfaceScalarField> zetasCorrOwn_;

        //- List of interpolated moments (neighbour)
        PtrList<surfaceScalarField> momentsNei_;

        //- List of interpolated moments (owner)
        PtrList<surfaceScalarField> momentsOwn_;

        //- Field to store the number of faces with outgoing flux per each cell
        mutable labelField nFacesOutgoingFlux_;

        //- Field to store the number of realizable moments in each cell
        mutable labelField nRealizableMoments_;

        //- Field to store the number of realizable m* in each cell
        mutable labelField nRealizableMomentsStar_;

        //- List of limiters for zeta_k
        PtrList<surfaceScalarField> limiters_;

        //- List of cell limiters
        PtrList<volScalarField> cellLimiters_;

        //- Face velocity
        const surfaceScalarField& phi_;


    // Private member functions

        //- Compute n values of zeta_k from n + 1 moments
        void computeZetaFields();

        //- Updates reconstructed moments from the corresponding values of zeta
        void updateMomentFieldsFromZetas
        (
            const surfaceScalarField& m0f,
            const PtrList<surfaceScalarField>& zetaf,
            PtrList<surfaceScalarField>& mf
        );

        //- Compute n + 1 moments from n values of zeta_k
        void zetaToMoments
        (
            const scalarList& zetaf,
            scalarList& mf,
            scalar m0 = 1.0
        );


protected:

    // Protected member functions

        //- Calculates the number of cells with outgoing flux
        void countFacesWithOutgoingFlux();

        //- Computes the limiter used for additional limitation
        void limiter();

        //- Reconstructs zeta_k
        void interpolateFields();

        //- Applies additional limitation to zeta_k, if needed
        void limitZetas();


public:

    //- Runtime type information
    TypeName("zeta");


    // Constructors

        //- Construct from univariateMomentSet
        zeta
        (
            const dictionary& dict,
            const scalarQuadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~zeta();


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
