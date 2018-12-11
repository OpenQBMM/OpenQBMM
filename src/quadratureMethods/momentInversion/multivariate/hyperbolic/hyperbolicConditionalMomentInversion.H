/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
     \\/     M anipulation  |
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
    Foam::hyperbolicConditionalMomentInversion

Description
    Implementation of the hyperbolic conditional quadrature method of moment
    for velocity. Three nodes are used in each direction. Moment realizability
    is verified analytically.

SourceFiles
    hyperbolicConditionalMomentInversion.C
    hyperbolicConditionalMomentInversionI.C

\*---------------------------------------------------------------------------*/

#ifndef hyperbolicConditionalMomentInversion_H
#define hyperbolicConditionalMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "vector2D.H"
#include "vector.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "mappedList.H"
#include "univariateMomentSet.H"
#include "multivariateMomentSet.H"
#include "hyperbolicMomentInversion.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class hyperbolicConditionalMomentInversion Declaration
\*---------------------------------------------------------------------------*/
class hyperbolicConditionalMomentInversion
{
    // Private member data

        //- Number of dimensions in velocity space
        const label nDimensions_;

        //- Number of moments
        const label nMoments_;

        //- List of the number of nodes
        const label nNodes_;

        //- Support of the NDF
        const word support_;

        //- Moments
        mappedList<scalar> moments_;

        //- Quadratrure abscissae
        mappedList<vector> abscissae_;

        //- Quadratrure weights
        mappedList<scalar> weights_;

        //- Univariate moment inversion method
        autoPtr<hyperbolicMomentInversion> univariateInverter_;

        //- Minimum value of the variance to compute eta
        scalar etaMin_;

        //- Maximum normalized skewness
        scalar qMax_;

        //- Small negative value for realizability
        scalar smallNegRealizability_;

        //- Minimum value of variance
        scalar varMin_;

        //- Minimum value of correlation
        scalar minCorrelation_;


    // Private member functions

        // Calculate q from q and eta
        scalar calcQ(scalar q, scalar eta);

        // Check realizability of univariate moments
        void realizabilityUnivariateMoments
        (
            scalar& c2,
            scalar& c3,
            scalar& c4
        );

        //- One-dimensional HyQMOM inversion
        void invert1D
        (
            const multivariateMomentSet& moments,
            scalarList& weights1D,
            scalarList& abscissae1D
        );

        //- Two-dimensional HyCQMOM inversion
        void invert2D
        (
            const multivariateMomentSet& moments,
            mappedList<scalar>& weights2D,
            mappedList<vector2D>& abscissae2D
        );

        //- Three-dimensional HyCQMOM inversion
        void invert3D
        (
            const multivariateMomentSet& moments
        );

        //- Reset inverter
        void reset();

        //- Return the number of moments given a number of dimensions
        label getNMoments(const label nDims) const;

        //- Return the moment orders give a number of dimensions
        labelListList getMomentOrders(const label nDims) const;

        //- Return the number of nodes given a number of dimensions
        label getNNodes(const label nDims) const;

        //- Return the node indexes give a number of dimensions
        labelListList getNodeIndexes(const label nDims) const;

        //- Disallow default bitwise copy construct
        //hyperbolicConditionalMomentInversion(const hyperbolicConditionalMomentInversion&);

        //- Disallow default bitwise assignment
        //void operator=(const hyperbolicConditionalMomentInversion&);


public:

    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        hyperbolicConditionalMomentInversion
        (
            const dictionary& dict,
            const label nDimensions
        );


    //- Destructor
    virtual ~hyperbolicConditionalMomentInversion();

    // Static data

        //- List of moment orders in two dimensions
        static const labelListList twoDimMomentOrders;

        //- List of moment orders in three dimensions
        static const labelListList threeDimMomentOrders;

        //- List of node indexes in two dimensions
        static const labelListList twoDimNodeIndexes;

        //- List of node indexes in three dimensions
        static const labelListList threeDimNodeIndexes;


    // Member Functions

        //- Invert moments to find weights and  abscissae
        void invert(const multivariateMomentSet& moments);

        //- Return number of moments
        inline label nMoments();

        //- Return number of number of nodes in each dimension
        inline label nNodes();

        //- Return mapped quadrature weights for each dimension
        inline const mappedList<scalar>& weights() const;

        //- Return mapped quadrature abscissae for each dimension
        inline const mappedList<vector>& abscissae() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hyperbolicConditionalMomentInversionI.H"

#endif


// ************************************************************************* //
