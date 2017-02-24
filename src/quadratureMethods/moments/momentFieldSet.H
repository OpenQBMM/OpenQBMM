/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2015-03-10 Alberto Passalacqua: Templated class on the type of moment and of
                                quadrature node.
2015-06-13 Alberto Passalacqua: Introduced autoPtr to the PtrList of nodes to
                                improve initialization of nodes.
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

#include "PtrList.H"
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
    public PtrList<momentType>
{
    // Private data

        //- Name of the momentFieldSet
        const word name_;

        //- Reference to the quadrature nodes
        const autoPtr<PtrList<nodeType>>& nodes_;

        //- The number of dimensions of the distribution function
        const label nDimensions_;

        //- The number of moments
        const label nMoments_;

        //- Map from the label representation of the moment name to its index
        //  within the List of moments
        Map<label> momentMap_;


    // Private Member Functions

        //- Given the desired dimension, convert to a label.
        //  e.g., If nDimensions_ = 3;
        //  0 -> 100
        //  1 -> 10
        //  2 -> 1
        //inline label dimToCmptOrder(label dim) const;

        //- Disallow default bitwise copy construct
        momentFieldSet(const momentFieldSet&);

        //- Disallow default bitwise assignment
        void operator=(const momentFieldSet&);


public:

    // Constructors

        //- Construct from dictionary and the quadrature nodes
        momentFieldSet
        (
            const word& distributionName,
            const dictionary& dict,
            const fvMesh& mesh,
            const autoPtr<PtrList<nodeType>>& nodes
        );

        //- Construct from size, quadrature nodes, dimensions and moment map
        //  The list of moments is initialised but moment fields are not
        //  automatically allocated.
        momentFieldSet
        (
            const word& distributionName,
            const label nMoments,
            const autoPtr<PtrList<nodeType>>& nodes,
            const label nDimensions,
            const Map<label>& momentMap
        );


    //- Destructor
    ~momentFieldSet();


    // Member Functions

        // Access

            //- Const access to the zeroth moment
            inline const momentType& operator()() const;

            //- Non-const access to the zeroth moment
            inline momentType& operator()();

            //- Const access to the first order moment in dimension a
            inline const momentType& operator()(label a) const;

            //- Non-const access to the first order moment in dimension a
            inline momentType& operator()(label a);

            //- Const access to the second order moment in dimension a+b
            inline const momentType& operator()(label a, label b) const;

            //- Non-const access to the second order moment in dimension a+b
            inline momentType& operator()(label a, label b);

            //- Const access to the third order moment in dimension a+b+c
            inline const momentType& operator()
            (
                label a,
                label b,
                label c
            ) const;

            //- Non-const access to the third order moment in dimension a+b+c
            inline momentType& operator()(label a, label b, label c);

            //- Returns const access to the moment map
            inline const Map<label>& momentMap() const;

            //- Return the dimensionality of the distribution
            inline label nDimensions() const;

            //- Recalculate the moments from the quadrature nodes
            void update();

            //- Return the map label of a moment
            label map(label a, label b, label c) const
            {
                return momentMap_
                [
                    100*a
                  + 10*b
                  + c
                ];
            };
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
