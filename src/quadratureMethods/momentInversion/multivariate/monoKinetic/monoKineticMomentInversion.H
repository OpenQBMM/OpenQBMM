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
    Foam::multivariateMomentInversions::monoKinetic

Description
    Implementation of the monokinetic velocity distribution moment set
    inversion

SourceFiles
    monoKineticmonoKineticMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef monoKineticMomentInversion_H
#define monoKineticMomentInversion_H


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
                        Class monoKinetic Declaration
\*---------------------------------------------------------------------------*/
class monoKinetic
:
    public multivariateMomentInversion
{
    // Private member data

        //- Number of size moments
        const label nSizeMoments_;

        //- Moment orders of conditional velocity distribution
        const labelListList velocityMomentOrders_;

        //- Univariate moment inversion method
        autoPtr<univariateMomentInversion> sizeInverter_;


    // Private member functions

        //- Calculate the number of pure size moments
        label calcNSizeMoments(const labelListList& momentOrders);


public:

    //- Runtime type information
    TypeName("monoKinetic");


    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        monoKinetic
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    //- Destructor
    virtual ~monoKinetic();


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
