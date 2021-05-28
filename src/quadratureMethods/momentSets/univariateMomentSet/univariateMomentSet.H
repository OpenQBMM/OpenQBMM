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
    Foam::univariateMomentSet

Description
    Stores a generic moment set with the corresponding list of orders.

SourceFiles
    univariateMomentSet.C
    univariateMomentSet.H

\*---------------------------------------------------------------------------*/

#ifndef univariateMomentSet_H
#define univariateMomentSet_H

#include "momentSet.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class univariateMomentSet Declaration
\*---------------------------------------------------------------------------*/

class univariateMomentSet : public momentSet
{
    // Private data

        //- Alpha coefficients of the recurrence relation
        scalarList alpha_;

        //- Beta coefficients of the recurrence relation
        scalarList beta_;

        //- Zeta values to check moment realizability
        scalarList zeta_;

        //- Canonical moments associated to the moment vector
        scalarList canonicalMoments_;

        //- Position of first negative zeta
        label negativeZeta_;

        //- Bool to indicate if moment set is degenerate, with only the
        //  zero-order moment realizable
        bool degenerate_;

        //- Bool to indicate if the full moment set is realizable
        bool fullyRealizable_;

        //- Bool to indicate if a subset of the full moment set is realizable
        bool subsetRealizable_;

        //- Bool to indicate if the moment set is on the boundary of the
        //  moment space
        bool onMomentSpaceBoundary_;

        //- Number of realizable moments
        label nRealizableMoments_;


     // Private member functions

        //- Check moment realizability
        void checkRealizability
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Checks canonical moments for measures with support over [0,1]
        void checkCanonicalMoments
        (
            const scalarList& zeta,
            const label nZeta
        );


protected:

        // Protected data

            //- Bool to indicate moment realizability was checked
            bool realizabilityChecked_;


public:

    // Constructors

        //- Construct from number of moments and initial value
        univariateMomentSet
        (
            const label nMoments,
            const word& support,
            const scalar initValue = 0,
            const label nFixedQuadraturePoints = 0
        );

        //- Construct from scalarList
        univariateMomentSet
        (
            const scalarList& m,
            const word& support,
            const label nFixedQuadraturePoints = 0
        );

    //- Destructor
    virtual ~univariateMomentSet();


    // Member Functions

        //- Returns alpha coefficients of the recurrence relation
        inline const scalarList& alphaRecurrence() const;

        //- Returns alpha coefficients of the recurrence relation
        inline scalarList& alphaRecurrence();

        //- Returns beta coefficients of the recurrence relation
        inline const scalarList& betaRecurrence() const;

        //- Returns beta coefficients of the recurrence relation
        inline scalarList& betaRecurrence();

        //- Returns zeta_k values for the moment set
        inline scalarList& zetas();

        //- Returns the canonical moments if support is [0, 1]
        inline scalarList& canonicalMoments();

        //- Returns degenerate status of the moment set
        inline bool isDegenerate();

        //- Return true if the entire moment set is realizable
        inline bool isFullyRealizable
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Return true if some moments are realizable
        inline bool isRealizable
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Return true if a subset of moments is realizable
        inline bool isSubsetRealizable
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        //- Return true if the moment set is on the boundary of the moment space
        inline bool isOnMomentSpaceBoundary
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        inline labelListList makeUnivariateMomentOrders
        (
            const label nMoments
        );

        //- Return the number of realizable moments
        inline label nRealizableMoments
        (
            bool fatalErrorOnFailedRealizabilityTest = true
        );

        inline label negativeZeta();

        //- Update moments from the corresponding quadrature
        virtual void update
        (
            const scalarList& weights,
            const scalarList& abscissae
        );

        //- Const access to element
        inline const scalar& operator[](const label index) const;

        //- Non-const access to element, reinitializing realizaibility check
        inline scalar& operator[](const label index);

        //- Const access to the moment of order momentOrder
        inline const scalar& operator()(const label momentOrder) const;

        //- Non-const access to the moment of order momentOrder
        inline scalar& operator()(const label momentOrder);

        //- Resizes moment set
        void resize(const label newSize);

        //- Set size of moment set
        void setSize(const label newSize);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
