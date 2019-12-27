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
    Foam::multivariateMomentInversion

Description
    General class for multivariate moment inversions

SourceFiles
    multivariateMomentInversion.C
    multivariateMomentInversionI.H

\*---------------------------------------------------------------------------*/

#ifndef multivariateMomentInversion_H
#define multivariateMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarList.H"
#include "vectorList.H"
#include "dictionary.H"
#include "mappedLists.H"
#include "Vandermonde.H"
#include "multivariateMomentSet.H"
#include "runTimeSelectionTables.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class multivariateMomentInversion Declaration
\*---------------------------------------------------------------------------*/
class multivariateMomentInversion
{
protected:
    // Protected member data

        //- Number of dimensions in distribution
        const label nDistributionDims_;

        //- Number of dimensions in velocity space
        const label nvelocityDimensions_;

        //- Moment orders of full distribution
        const labelListList& momentOrders_;

        //- Node indexes of full distribution
        const labelListList& nodeIndexes_;

        //- List of velocity indexes
        labelList velocityIndexes_;

        //- Number of nodes in each direction
        labelList nNodes_;

        //- Quadratrure weights
        mappedScalarList weights_;

        //- Quadratrure size abscissae
        mappedList<scalarList> abscissae_;

        //- Quadratrure velocity abscissae
        mappedVectorList velocityAbscissae_;


    // Protected member functions

        //- Reset inverter
        void reset();

        //- Compare indexes up to the minimum size
        static bool compare(const labelList& index1, const labelList& index2);

        //- Compare indexes up to dimi
        static bool compare
        (
            const labelList& index1,
            const labelList& index2,
            const label size
        );


public:

    //- Runtime type information
    TypeName("multivariateMomentInversion");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        multivariateMomentInversion,
        dictionary,
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        ),
        (dict, momentOrders, nodeIndexes, velocityIndexes)
    );


    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        multivariateMomentInversion
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    // Selectors

        static autoPtr<multivariateMomentInversion> New
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    //- Destructor
    virtual ~multivariateMomentInversion();


    // Member Functions

        //- Invert moments to find weights and  abscissae
        virtual bool invert(const multivariateMomentSet& moments) = 0;

        //- Return mapped quadrature weights for each dimension
        inline const mappedScalarList& weights() const;

        //- Return mapped quadrature size abscissae for each dimension
        inline const mappedList<scalarList>& abscissae() const;

        //- Return mapped quadrature velocity abscissae for each dimension
        inline const mappedVectorList& velocityAbscissae() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multivariateMomentInversionI.H"

#endif


// ************************************************************************* //
