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
    Foam::multivariateMomentInversions::TensorProduct

Description
    Tensor product moment inversion for multidimensional moment sets

SourceFiles
    TensorProductMomentInversion.C

\*---------------------------------------------------------------------------*/

#ifndef TensorProductMomentInversion_H
#define TensorProductMomentInversion_H


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
                        Class TensorProduct Declaration
\*---------------------------------------------------------------------------*/
class TensorProduct
:
    public multivariateMomentInversion
{
    // Private member data

        //- Number of pure moments in each direction
        labelList nPureMoments_;

        //- Support of the NDF
        wordList supports_;

        //- Univariate moment inversion method
        PtrList<univariateMomentInversion> univariateInverters_;


    // Private member functions

        //- Construct non-zero nodes indexes
        void buildIndexes
        (
            labelListList& nodeIndexes,
            const labelList& nNodes,
            label dimi,
            label& nodei,
            labelList& index
        );


public:

    //- Runtime type information
    TypeName("TensorProduct");


    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        TensorProduct
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes,
            const labelList& velocityIndexes
        );


    //- Destructor
    virtual ~TensorProduct();


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
