/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2012-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 Alberto Passalacqua
-------------------------------------------------------------------------------
2015-03-10 Alberto Passalacqua: Templated class on the type of moment and of
                                quadrature node.
2015-06-13 Alberto Passalacqua: Introduced autoPtr to the PtrList of nodes to
                                improve initialization of nodes.
2017-03-26 Alberto Passalacqua: Added the capability to recompute moments
                                locally.
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
    Foam::momentFieldSet

Description
    Stores a list of moments, with a Map<label> that takes a label key and
    returns the index into the list of moments.

    Example

        Moment  |  Label Key | Index
         0 0 0         0         0
         1 0 0       100         1
         0 1 0        10         2
         0 0 1         1         3
         2 0 0       200         4
         2 1 0       210         5

    operator() is overloaded in order to hide the above mapping.

SourceFiles
    momentFieldSetI.H
    momentFieldSet.C

\*---------------------------------------------------------------------------*/

#ifndef momentFieldSet_H
#define momentFieldSet_H

#include "mappedPtrList.H"
#include "volFields.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class momentFieldSet Declaration
\*---------------------------------------------------------------------------*/

template <class momentType, class nodeType>
class momentFieldSet
:
    public mappedPtrList<momentType>
{
    // Private data

        //- Name of the momentFieldSet
        const word name_;

        //- Reference to the quadrature nodes
        const autoPtr<mappedPtrList<nodeType>>& nodes_;

        //- The number of dimensions of the distribution function
        const label nDimensions_;

        //- The number of moments
        const label nMoments_;

        //- Support of the measure associated to the moment moment set
        word support_;


    // Private Member Functions

        //- Given the desired dimension, convert to a label.
        //  e.g., If nDimensions_ = 3;
        //  0 -> 100
        //  1 -> 10
        //  2 -> 1
        //inline label dimToCmptOrder(label dim) const


public:

    // Constructors

        //- Construct from dictionary and the quadrature nodes
        momentFieldSet
        (
            const word& distributionName,
            const dictionary& dict,
            const fvMesh& mesh,
            const autoPtr<mappedPtrList<nodeType>>& nodes,
            const word& support
        );

        //- Construct from size, quadrature nodes, dimensions and moment map
        //  The list of moments is initialized but moment fields are not
        //  automatically allocated.
        momentFieldSet
        (
            const word& distributionName,
            const label nMoments,
            const autoPtr<mappedPtrList<nodeType>>& nodes,
            const label nDimensions,
            const Map<label>& momentMap,
            const word& support
        );


        //- Disallow default bitwise copy construct
        momentFieldSet(const momentFieldSet&) = delete;


    //- Destructor
    ~momentFieldSet();


    // Member Functions

        // Access

            //- Return the dimensionality of the distribution
            inline label nDimensions() const;

            //- Return support of the measure associated to the moment set
            inline const word support() const;

            //- Recalculate the moments from the quadrature nodes
            void update();

            //- Recalculate the moments from the quadrature nodes
            //  on boundaries
            void updateBoundaries();

            //- Recalculate the moments from the quadrature nodes in an element
            void updateLocalMoments(label elemi);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const momentFieldSet&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "momentFieldSetI.H"
#   include "momentFieldSet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
