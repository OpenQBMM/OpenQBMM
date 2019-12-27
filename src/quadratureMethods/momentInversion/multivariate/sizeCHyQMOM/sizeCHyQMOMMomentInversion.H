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
    Foam::multivariateMomentInversions::sizeCHyQMOM

Description
    Implementation of the hyperbolic conditional quadrature method of moment
    for velocity with conditioning on size.

SourceFiles
    sizeCHyQMOMMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef sizeCHyQMOMMomentInversion_H
#define sizeCHyQMOMMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "vector2D.H"
#include "vector.H"
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
                        Class sizeCHyQMOM Declaration
\*---------------------------------------------------------------------------*/
template<class velocityInversion>
class sizeCHyQMOMBase
:
    public multivariateMomentInversion
{
    // Private member data

        //- Number of pure size moments
        const label nSizeMoments_;

        //- Moment orders of conditional velocity distribution
        const labelListList velocityMomentOrders_;

        //-Number of nodes in size direction
        const label nSizeNodes_;

        //- Node indexes of conditional velocity distribution
        const labelListList velocityNodeIndexes_;

        //- Univariate moment inversion method
        autoPtr<univariateMomentInversion> sizeInverter_;

        //- Velocity moment inversion method
        autoPtr<velocityInversion> velocityInverter_;


    // Private member functions

        //- Calculate the number of pure size moments
        label calcNSizeMoments(const labelListList& momentOrders);


public:

    //- Runtime type information
    TypeName("sizeCHyQMOMBase");

    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        sizeCHyQMOMBase
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    //- Destructor
    virtual ~sizeCHyQMOMBase();


    // Member Functions

        //- Invert moments to find weights and  abscissae
        virtual bool invert(const multivariateMomentSet& moments);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace multivariateMomentInversions
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sizeCHyQMOMMomentInversion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
