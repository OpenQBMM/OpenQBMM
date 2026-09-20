/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2026 Alberto Passalacqua
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
    Foam::multivariateMomentInversions::conditionalMeans

Description
    Inversion of a distribution in one primary coordinate and any number of
    coordinates conditioned on it, each carried by its mean at every node
    of the primary coordinate.

    The moments of the primary coordinate, m_k = <x^k> for k < 2n, are
    inverted to n nodes with a univariate inversion, and the conditional
    mean of every other coordinate phi at each node is found from the n
    moments <x^k phi>, k < n, through the shared conditional Vandermonde
    solve. Nothing of the conditioned coordinates beyond their means is
    carried: the distribution is a delta in them at every node.

    The primary coordinate is the first dimension of the distribution that
    is not a velocity. A conditioned coordinate is a velocity, stored with
    the velocity abscissae as a monokinetic closure does, or a scalar,
    stored with the abscissae, as the compositions conditioned on the
    mixture fraction of a mixing-reaction closure are.

    The support of every coordinate comes from the moment set, which takes
    it from the model the quadrature belongs to. The primary moments are
    inverted, and tested for realizability, on theirs. A conditioned scalar
    on the positive reals or on [0, 1] is clipped to it when the solve
    gives a mean just outside, which transport error can do; how many
    values were clipped in the last inversion is kept, so that a caller can
    report it. A velocity is on the reals and is never clipped.

    Clipping a value is a repair of one node on its own, which changes the
    moment it belongs to: harmless where the value is outside its support
    by round-off, and not where it is by more. The conditional means are
    what reproduces the moments, which is the value of the scalar at the
    nodes only where the scalar is a polynomial of the primary coordinate
    of a degree below the number of nodes; a profile with a corner, which a
    fast reaction gives a composition, comes out below zero on one side of
    it. A model that has a better repair, as a reaction has along its
    stoichiometry, asks the inversion not to clip, and receives the
    conditional means as the solve gives them.

    Selected as conditionalMeans, or as monoKinetic, the name the velocity
    closure has always had.

SourceFiles
    conditionalMeansMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef conditionalMeansMomentInversion_H
#define conditionalMeansMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "vectorList.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "mappedLists.H"
#include "Vandermonde.H"
#include "multivariateMomentInversion.H"
#include "univariateMomentInversion.H"


namespace Foam
{
namespace multivariateMomentInversions
{

/*---------------------------------------------------------------------------*\
                      Class conditionalMeans Declaration
\*---------------------------------------------------------------------------*/

class conditionalMeans
:
    public multivariateMomentInversion
{
    // Private member data

        //- Index of the primary dimension in the moment orders: the first
        //  dimension of the distribution that is not a velocity. It is not
        //  required to be the first dimension: the size of a monokinetic
        //  closure may come after its velocities.
        const label primaryIndex_;

        //- The dimensions of the conditioned scalars, and the component of
        //  the abscissae each is stored in
        labelList scalarIndexes_;
        labelList scalarComponents_;

        //- Number of moments of the primary coordinate
        const label nPrimaryMoments_;

        //- Univariate inversion of the primary coordinate
        autoPtr<univariateMomentInversion> primaryInverter_;

        //- Number of conditioned values clipped to their support in the
        //  last inversion
        label nClipped_;

        //- Whether a conditioned scalar is kept within its support
        bool clipToSupports_;


    // Private member functions

        //- The index of the primary dimension, stopping if every dimension
        //  is a velocity
        static label calcPrimaryIndex
        (
            const label nDistributionDims,
            const labelList& velocityIndexes
        );

        //- The number of pure moments of the primary coordinate
        static label calcNPrimaryMoments
        (
            const labelListList& momentOrders,
            const label primaryIndex
        );


public:

    //- Runtime type information
    TypeName("conditionalMeans");


    // Constructors

        //- Construct from the dictionary of the inversion, the moment
        //  orders, the node indexes and the velocity dimensions
        conditionalMeans
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    //- Destructor
    virtual ~conditionalMeans();


    // Static Member Functions

        //- The moment orders and node indexes of a distribution with
        //  nPrimaryMoments moments of its primary coordinate, conditioned
        //  on which are nScalars scalars and then nVelocityDimensions
        //  velocity components, in this order of dimensions: the pure
        //  moments of the primary coordinate, then the first-order moment
        //  of each conditioned coordinate times the powers of the primary
        //  coordinate below the number of its nodes.
        static void momentOrdersAndNodes
        (
            const label nPrimaryMoments,
            const label nScalars,
            const label nVelocityDimensions,
            labelListList& momentOrders,
            labelListList& nodeIndexes
        );


    // Member Functions

        //- Invert moments to find weights and abscissae
        virtual bool invert(const multivariateMomentSet& moments);

        //- Whether a conditioned scalar is kept within its support
        virtual void clipToSupports(const bool clip)
        {
            clipToSupports_ = clip;
        }

        //- Number of conditioned values clipped to their support in the
        //  last inversion
        virtual label nClipped() const
        {
            return nClipped_;
        }

        //- Return the smallest accepted value for m0
        virtual scalar smallM0() const;

        //- Return the smallest accepted value for zeta_k
        virtual scalar smallZeta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multivariateMomentInversions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
