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
    // Private member functions

        //- Disallow default bitwise copy construct
        univariateMomentInversion(const univariateMomentInversion&);

        //- Disallow default bitwise assignment
        void operator=(const univariateMomentInversion&);


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
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from univariateMomentSet
        univariateMomentInversion
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~univariateMomentInversion();


    // Selectors

        static autoPtr<univariateMomentInversion> New
        (
            const dictionary& dict
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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateMomentInversionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //