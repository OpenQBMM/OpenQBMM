/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
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
    Foam::multivariateMomentInversions::conditional

Description
    Implementation of the conditional quadrature method of moments

    References
    \verbatim
        "Conditional quadrature method of moments for kinetic equations"
        C Yuan, R O Fox
        Journal of Computational Physics
        Volume 230, Pages 8216-8246, 2011
    \endverbatim

SourceFiles
    conditionalMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef conditionalMomentInversion_H
#define conditionalMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "multivariateMomentInversion.H"
#include "univariateMomentSet.H"
#include "univariateMomentInversion.H"
#include "multivariateMomentSet.H"
#include "Vandermonde.H"
#include "mappedLists.H"
#include "PtrList.H"

namespace Foam
{
namespace multivariateMomentInversions
{

/*---------------------------------------------------------------------------*\
                     Class conditionalMomentInversion Declaration
\*---------------------------------------------------------------------------*/
class conditional
:
    public multivariateMomentInversion
{
    // Private member data

        //- Number of pure moments in each direction
        labelList nPureMoments_;

        //- Type of support for each direction
        wordList supports_;

        //- Stored moments
        multivariateMomentSet moments_;

        //- Conditional weights
        PtrList<mappedList<scalar>> conditionalWeights_;

        //- Stored conditional moments
        PtrList<PtrList<mappedList<scalar>>> conditionalMoments_;

        //- Inverted VanderMonde matricies
        PtrList<mappedList<scalarSquareMatrix> > invVR_;

        //- Univariate moment inversion method
        PtrList<univariateMomentInversion> momentInverters_;

        //- Current indexes
        label vi_;
        label si_;


    // Private member functions

        //- Create map for nodes
        void setNodeMap
        (
            Map<label>& map,
            const label nDimensions,
            const labelList& nNodes,
            label dimi,
            label& mi,
            labelList& pos
        );

        //- Cycles through all alphas and solves the given conditional moments
        void cycleAlphaCM
        (
            const label dimi,
            const label dimJ,
            label ai,
            labelList& pos
        );

        //- Sets the Vandermonde matricies for the previous dimension
        void setVR
        (
            const label dimi,
            labelList& pos,
            label ai
        );

        //- Once all of the conditional moments are known, said moments are
        //  inverted to find the weights and abscissae
        bool cycleAlphaWheeler
        (
            const label dimi,
            label alphai,
            labelList& pos
        );


public:

    //- Runtime type information
    TypeName("conditional");


    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        conditional
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    //- Destructor
    virtual ~conditional();


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
