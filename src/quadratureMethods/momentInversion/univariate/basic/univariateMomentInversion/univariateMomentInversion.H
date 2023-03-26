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
    Copyright (C) 2019-2023 Alberto Passalacqua
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
    Foam::univariateMomentInversion

Description
    Abstract class for univariate quadrature on a moment set.

SourceFiles
    univariateMomentInversion.C
    univariateMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef univariateMomentInversion_H
#define univariateMomentInversion_H

#include "scalarList.H"
#include "scalarMatrices.H"
#include "runTimeSelectionTables.H"
#include "univariateMomentSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class univariateMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class univariateMomentInversion
{
    // Private data

        //- Smallest accepted value for m0
        scalar smallM0_;

        //- Smallest accepted value for zeta_k
        scalar smallZeta_;


protected:

    // Protected data

        //- Number of invertible moments for the user-selected quadrature
        label nInvertibleMoments_;

        //- Number of quadrature nodes
        label nNodes_;

        //- Quadrature abscissae
        scalarList abscissae_;

        //- Quadrature weights
        scalarList weights_;


    // Protected member functions

        //- Create Jacobi matrix
        void JacobiMatrix
        (
            univariateMomentSet& moments,
            scalarSquareMatrix& z,
            const scalar minKnownAbscissa = 0,
            const scalar maxKnownAbscissa = 0
        );


        // Calculate weights and abscissae
        void calcQuadrature
        (
            const univariateMomentSet& moments,
            const scalarSquareMatrix& z
        );


public:

    //- Runtime type information
    TypeName("univariateMomentInversion");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            univariateMomentInversion,
            dictionary,
            (
                const dictionary& dict,
                const label nMaxNodes
            ),
            (dict, nMaxNodes)
        );


    // Constructors

        //- Construct from univariateMomentSet
        univariateMomentInversion
        (
            const dictionary& dict,
            const label nMaxNodes = 0
        );

        //- Disallow default bitwise copy construct
        univariateMomentInversion(const univariateMomentInversion&) = delete;


    //- Destructor
    virtual ~univariateMomentInversion();


    // Selectors

        static autoPtr<univariateMomentInversion> New
        (
            const dictionary& dict,
            const label nMaxNodes = 0
        );


    // Member Functions

        //- Calculates the number of quadrature nodes
        virtual void calcNQuadratureNodes
        (
            univariateMomentSet& moments
        ) = 0;

        //- Modify the recurrence relation for fixed-point quadrature
        virtual void correctRecurrence
        (
            univariateMomentSet& moments,
            scalarList& alpha,
            scalarList& beta,
            const scalar minKnownAbscissa = 0,
            const scalar maxKnownAbscissa = 0
        ) = 0;

        //- Invert the set of moments to compute weights and abscissae
        virtual void invert
        (
            univariateMomentSet& moments,
            const scalar minKnownAbscissa = 0,
            const scalar maxKnownAbscissa = 0
        );

        //- Return quadrature abscissae
        inline const scalarList& abscissae() const;

        //- Return the number of nodes
        inline label nNodes() const;

        //- Return quadrature weigths
        inline const scalarList& weights() const;

        //- Return the smallest accepted value for m0
        inline scalar smallM0() const;

        //- Return the smallest accepted value for zeta_k
        inline scalar smallZeta() const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const univariateMomentInversion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentInversionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
