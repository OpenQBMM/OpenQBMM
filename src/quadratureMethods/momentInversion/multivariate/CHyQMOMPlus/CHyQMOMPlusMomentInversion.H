/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 Alberto Passalacqua
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
    Foam::multivariateMomentInversions::CHyQMOMPlus

Description
    Implementation of the hyperbolic conditional quadrature method of moment
    for velocity. Three nodes are used in each direction. Moment realizability
    is verified analytically.

SourceFiles
    hyperbolicConditionalMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef CHyQMOMPlusMomentInversion_H
#define CHyQMOMPlusMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "vector2D.H"
#include "vector.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "mappedList.H"
#include "univariateMomentSet.H"
#include "multivariateMomentInversion.H"
#include "hyperbolicMomentInversion.H"


namespace Foam
{
namespace multivariateMomentInversions
{

/*---------------------------------------------------------------------------*\
                    Class CHyQMOMPlus Declaration
\*---------------------------------------------------------------------------*/
class CHyQMOMPlus
:
    public multivariateMomentInversion
{
    // Private member data

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

        //- Disallow default bitwise copy construct
        //CHyQMOMPlus(const CHyQMOMPlus&);

        //- Disallow default bitwise assignment
        //void operator=(const CHyQMOMPlus&);


public:

    //- Runtime type information
    TypeName("CHyQMOMPlus");


    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        CHyQMOMPlus
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    //- Destructor
    virtual ~CHyQMOMPlus();

    // Static data

        //- Return the number of moments given a number of dimensions
        static label getNMoments(const label nDims);

        //- Return the moment orders give a number of dimensions
        static labelListList getMomentOrders(const label nDims);

        //- Return the number of nodes given a number of dimensions
        static label getNNodes(const label nDims);

        //- Return the node indexes give a number of dimensions
        static labelListList getNodeIndexes(const label nDims);

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
        virtual bool invert(const multivariateMomentSet& moments);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multivariateMomentInversions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
