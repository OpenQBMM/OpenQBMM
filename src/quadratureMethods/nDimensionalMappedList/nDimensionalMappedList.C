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

\*---------------------------------------------------------------------------*/

#include "nDimensionalMappedList.H"


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template <class mappedType>
Foam::word
Foam::nDimensionalMappedList<mappedType>::listToWord(const labelList& lst)
{
    word w;

    forAll(lst, dimi)
    {
        w += Foam::name(lst[dimi]);
    }

    return w;
}

template <class mappedType>
Foam::label
Foam::nDimensionalMappedList<mappedType>::listToLabel(const labelList& lst)
{
    label l = 0;

    forAll(lst, dimi)
    {
        l += lst[dimi]*pow(10, lst.size() - dimi - 1);
    }

    return l;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class mappedType>
Foam::nDimensionalMappedList<mappedType>::nDimensionalMappedList
(
    const label nCmpt,
    const label nDims,
    const Map<label>& map
)
:
    PtrList<mappedType>(nCmpt),
    nDims_(nDims),
    map_(map)
{}

template <class mappedType>
Foam::nDimensionalMappedList<mappedType>::nDimensionalMappedList
(
    const label nDims,
    const labelList& nNodes
)
:
    PtrList<mappedType>(nDimensionalListLength(nDims, nNodes)),
    nDims_(nDims),
    map_(this->size())
{
    labelList pos(nDims);
    label mi = 0;
    setMappedPositions(nNodes, 0, mi, pos);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class mappedType>
Foam::nDimensionalMappedList<mappedType>::~nDimensionalMappedList()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class mappedType>
void Foam::nDimensionalMappedList<mappedType>::setMappedPositions
(
    const labelList& nNodes,
    label dimi,
    label& mi,
    labelList& pos
)
{
    if (dimi < nDims_)
    {
        for (label i = 0; i < nNodes[dimi]; i++)
        {
            pos[dimi] = i;
            setMappedPositions(nNodes, dimi+1, mi, pos);
        }
    }
    else
    {
        map_.insert
        (
            listToLabel
            (
                pos
            ),
            mi
        );
        mi++;
    }
}

template <class mappedType>
Foam::label Foam::nDimensionalMappedList<mappedType>::nDimensionalListLength
(
    const label nDims,
    const labelList nNodes
)
{
    label totalProduct = 1;
    for(label nodei = 0; nodei < nDims; nodei++)
    {
        totalProduct *= nNodes[nodei];
    }
    return totalProduct;
}


// ************************************************************************* //
