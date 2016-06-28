/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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
    Foam::nDimensionalMappedList

Description
    Stores a list, with a Map<label> that takes a label key and
    returns the index into the list.

    Example

        Lookup  |  Label Key | Index
         0 0 0         0         0
         1 0 0       100         1
         0 1 0        10         2
         0 0 1         1         3
         2 0 0       200         4
         2 1 0       210         5

    operator() is overloaded in order to hide the above mapping.

SourceFiles
    nDimensionalMappedListI.H
    nDimensionalMappedList.C

\*---------------------------------------------------------------------------*/

#ifndef nDimensionalMappedList_H
#define nDimensionalMappedList_H

#include "scalar.H"
#include "PtrList.H"
#include "volFields.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class nDimensionalMappedList Declaration
\*---------------------------------------------------------------------------*/

template <class mappedType>
class nDimensionalMappedList
:
    public PtrList<mappedType>
{
    // Private data

        //- The number of dimensions for the mapped data set
        const label nDims_;

        //- Map from the label representation of the moment name to its index
        //  within the List
        Map<label> map_;


    // Private Member Functions

        //- Given the desired dimension, and nNodes, create a map.
        //  e.g., If nNodes = [2,2];
        //  0 -> 00
        //  1 -> 10
        //  2 -> 01
        //  3 -> 11
        inline void setMappedPositions
        (
            const labelList& nNodes,
            label dimI,
            label& mI,
            labelList& pos
        );

        //- Calculate the length of the PtrList given nNodes, and nDimensions
        label nDimensionalListLength
        (
            const label nDims,
            const labelList nNodes
        );

        //- Disallow default bitwise copy construct
        nDimensionalMappedList(const nDimensionalMappedList&);

        //- Disallow default bitwise assignment
        void operator=(const nDimensionalMappedList&);


public:

    // Constructors

        //- Construct from size, dimensions, and map
        //  The list is initialised
        nDimensionalMappedList
        (
            const label nCmpt,
            const label nDimensions,
            const Map<label>& momentMap
        );

        //- Construct mapped list using the number of nodes in each direction,
        //  nDims is not necessarily equal to nNodes.size()
        nDimensionalMappedList
        (
            const label nDims,
            const labelList& nNodes
        );

    //- Destructor
    ~nDimensionalMappedList();


    // Member Functions

        // Static member functions
            //- Convert a list of labels to a word. {1, 2, 3} -> 123
            static word listToWord(const labelList& lst);

            //- Convert a list of labels to a single label. {0, 2, 3} -> 23
            static label listToLabel(const labelList& lst);


        // Access

            //- Const access to the first element
            inline const mappedType& operator()() const;

            //- Non-const access to the first element
            inline mappedType& operator()();

            //- Const access to the first dimension
            inline const mappedType& operator()(label a) const;

            //- Non-const access to the first dimension
            inline mappedType& operator()(label a);

            //- Const access to the second dimension
            inline const mappedType& operator()(label a, label b) const;

            //- Non-const access to the second dimension
            inline mappedType& operator()(label a, label b);

            //- Const access to the third dimension
            inline const mappedType& operator()
            (
                label a,
                label b,
                label c
            ) const;

            //- Non-const access to the third dimension
            inline mappedType& operator()
            (
                label a,
                label b,
                label c
            );

            //- Const access to the fourth dimension
            inline const mappedType& operator()
            (
                label a,
                label b,
                label c,
                label d
            ) const;

            //- Non-const access to the fourth dimension
            inline mappedType& operator()
            (
                label a,
                label b,
                label c,
                label d
            );

            //- Const access to the fith dimension
            inline const mappedType& operator()
            (
                label a,
                label b,
                label c,
                label d,
                label e
            ) const;

            //- Non-const access to the fith dimension
            inline mappedType& operator()
            (
                label a,
                label b,
                label c,
                label d,
                label e
            );

            //- Constant access given a label list
            inline const mappedType& operator()(const labelList& l) const;

            //- Non-constant access given a label list
            inline mappedType& operator()(const labelList& l);

            //- Returns const access to the moment map
            inline const Map<label>& map() const;

            //- Return the number of dimensions of the data set
            inline label nDims() const;

            //- Return the map label of a lookup
            label map(label a, label b, label c, label d, label e) const
            {
                return map_
                [
                    10000*a
                  + 1000*b
                  + 100*c
                  + 10*d
                  + e
                ];
            };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "nDimensionalMappedListI.H"
#   include "nDimensionalMappedList.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
