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
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    Second-order univariate moment advection with zeta simplified scheme for
    moments of NDF supported on R+ or [0, 1].

    Auxiliary quantities are used to perform the reconstruction. These 
    quantities are the zeta_k scalars determined when verifying the 
    realizability of the moment vector, in the case of support over R+. 
    Canonical moments are used when suppor is over [0, 1].

    This scheme is suitable, for example, for population balance and turbulent 
    mixing applications.

    References
    \verbatim
        "Realizable second-order finite-volume schemes for the advection of
        moment sets of the particle size distribution"
        F Laurent, T. T. Nguyen
        Journal of Computational Physics
        Volume 337, Pages 309-338, 2017
    \endverbatim

    \verbatim
        "A second-order realizable scheme for moment advection on unstructured 
        grids"
        A Passalacqua, F Laurent, R. O. Fox
        Computer Physics Communications
        Volume 248, Page 106993, 2020
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

        //- Number of zeta_k or p_k values
        label nAuxiliaryFields_;

        //- List of fields of zeta_k or p_k (n fields for n + 1 moments)
        PtrList<volScalarField> auxiliaryFields_;

        //- List of interpolated zeta_k or p_k values (neighbour)
        PtrList<surfaceScalarField> auxiliaryFieldsNei_;

        //- List of interpolated nodes (owner)
        PtrList<surfaceScalarField> auxiliaryFieldsOwn_;

        //- List of interpolated zeta_k or p_k values (neighbour)
        PtrList<surfaceScalarField> auxiliaryFieldsUpwindNei_;

        //- List of interpolated zeta_k or p_k (owner)
        PtrList<surfaceScalarField> auxiliaryFieldsUpwindOwn_;

        //- List of second-order correction for zeta_k or p_k values (neighbour)
        PtrList<surfaceScalarField> auxiliaryFieldsCorrNei_;

        //- List of second-order correction for zeta_k or p_k (owner)
        PtrList<surfaceScalarField> auxiliaryFieldsCorrOwn_;

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

        //- List of limiters for auxiliary fields
        PtrList<surfaceScalarField> limiters_;

        //- List of cell limiters
        PtrList<volScalarField> cellLimiters_;

        //- Face velocity
        const surfaceScalarField& phi_;


    // Private member functions

        //- Compute n values of auxiliary fields from n + 1 moment fields
        void computeAuxiliaryFields();

        //- Updates reconstructed moments from the corresponding values of
        //  auxiliary quantities
        void updateMomentFieldsFromAuxiliaryQuantities
        (
            const surfaceScalarField& m0f,
            const PtrList<surfaceScalarField>& auxiliaryFieldsf,
            PtrList<surfaceScalarField>& mf
        );

        //- Compute n + 1 moments from n values of zeta_k
        void zetaToMoments
        (
            const scalarList& zetaf,
            scalarList& mf,
            scalar m0 = 1.0
        );

        //- Compute n + 1 moments from n canonical moments
        void canonicalMomentsToMoments
        (
            const scalarList& canononicalMomentsf,
            scalarList& mf,
            scalar m0
        );

        //- Compute n + 1 moments from n values of the auxiliary quantity
        void auxiliaryQuantitiesToMoments
        (
            const scalarList& auxiliaryQuantityf,
            scalarList& mf,
            scalar m0 = 1.0
        );

protected:

    // Protected member functions

        //- Calculates the number of cells with outgoing flux
        void countFacesWithOutgoingFlux();

        //- Reconstructs auxiliary fields
        void interpolateFields();

        //- Applies additional limitation to auxiliary fields, if needed
        void limitAuxiliaryFields();


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
