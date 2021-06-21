/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    Foam::mappedPtrList

Description
    Stores a pointer list, with a Map<label> that takes a label key and
    returns the index into the list.

    Example

        Lookup  |  Label Key | Index
         0 0 0         0         0
         1 0 0       100         1
         0 1 0        10         2
         0 0 1         1         3
         2 0 0       200         4
         2 1 0       210         5

    The operator () is overloaded in order to hide the above mapping.

SourceFiles
    mappedPtrListI.H
    mappedPtrList.C

\*---------------------------------------------------------------------------*/

#ifndef mappedPtrList_H
#define mappedPtrList_H

#include "scalar.H"
#include "PtrList.H"
#include "labelList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mappedList Declaration
\*---------------------------------------------------------------------------*/

template <class mappedType>
class mappedPtrList
:
    public PtrList<mappedType>
{
    // Private data

        //- Map from the label to its index within the List
        Map<label> map_;

        //- Number of dimensions
        label nDimensions_;

    // Private Member Functions

        //- Returns an element provided the indeces of the order
        label calcMapIndex(std::initializer_list<label> indexes) const;


public:

    // Constructors

        //- Construct from size and labelListList
        mappedPtrList
        (
            const label size,
            const labelListList& indexes
        );

        //- Construct from size and map with initial value for all elements
        mappedPtrList
        (
            const label size,
            const Map<label>& map
        );

        //- Construct from List<mappedType> and labelListList
        mappedPtrList
        (
            const PtrList<mappedType>& initList,
            const labelListList& indexes
        );

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        mappedPtrList(Istream&, const INew&);


    //- Destructor
    ~mappedPtrList();


    // Member Functions

        // Static member functions

            //- Convert a list of labels to a word. {1, 2, 3} -> 123
            static word listToWord(const labelList& list);

            //- Convert a list of labels to a single label. {0, 2, 3} -> 23
            static label listToLabel
            (
                const labelList& list,
                const label nDimensions = 0
            );


        // Access functions

            //- Non-const access to a moment given the list of orders
            template <typename ...ArgsT>
            inline mappedType& operator()(ArgsT ...args);

            //- Const access to a moment given the list of orders
            template <typename ...ArgsT>
            inline const mappedType& operator()(ArgsT ...args) const;

            //- Constant access given a label list
            inline const mappedType& operator()(const labelList& list) const;

            //- Non-constant access given a label list
            inline mappedType& operator()(const labelList& list);

            //- Returns const access to the moment map
            inline const Map<label>& map() const;


        // Return functions

            //- Is index i set
            bool set(const label i) const;

            //- Is label list index set
            bool set(const labelList& list) const;

            //- Is label list index used
            bool found(const labelList& list) const;

            //- Is label list index used
            template <typename ...ArgsT>
            bool found(ArgsT ...args) const;


        // Edit functions

            //- Set map after PtrList has been constructed
            void setMap(const Map<label>& map);

            //- Set mapped location using pointer
            void set(const label i, mappedType* entry);

            //- Set mapped location using pointer and label list
            void set(const labelList& list, mappedType* entry);

            //- Set mapped location using autoPtr and label list
            void set(const labelList& list, autoPtr<mappedType> entry);

            //- Set mapped location using tmp and label list
            void set(const labelList& list, tmp<mappedType> entry);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mappedPtrListI.H"
#   include "mappedPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
