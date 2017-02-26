/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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
#include "univariateMomentSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class univariateMomentInversion Declaration
\*---------------------------------------------------------------------------*/

class univariateMomentInversion
{
protected:

    // Protected data

        //- Set of moments to invert
        univariateMomentSet& moments_;

        //- Bool to indicate if moments were inverted
        bool inverted_;

        //- Number of invertible moments
        label nInvertibleMoments_;

        //- Number of quadrature nodes
        label nNodes_;

        //- Quadrature weights
        scalarList weights_;

        //- Quadrature abscissae
        scalarList abscissae_;


    // Protected member functions

        //- Returns the list of alpha parameters for the quadrature
        virtual scalarList& alphaRecurrence();

        //- Returns the list of beta parameters for the quadrature
        virtual scalarList& betaRecurrence();

        //- Calculates the number of quadrature nodes
        virtual void calcNQuadratureNodes() = 0;

        //- Create Jacobi matrix
        scalarSquareMatrix JacobiMatrix
        (
            const scalarList& alpha,
            const scalarList& beta
        );


public:

    // Constructors

        //- Construct from univariateMomentSet
        univariateMomentInversion
        (
            univariateMomentSet& moments
        );

    //- Destructor
    virtual ~univariateMomentInversion();


    // Member Functions

        //- Return quadrature abscissae
        inline const scalarList& abscissae() const;

        //- Invert the set of moments to compute weights and abscissae
        virtual void invert();

        //- Returns the status of inversion
        inline bool isInverted() const;

        //- Return the number of invertible moments
        inline label nInvertibleMoments() const;

        //- Return the number of nodes
        inline label nNodes() const;

        //- Return quadrature weights
        inline const scalarList& weights() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentInversionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
