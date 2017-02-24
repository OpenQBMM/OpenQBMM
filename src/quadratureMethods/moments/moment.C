/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2015-03-09 Alberto Passalacqua: Templated class on the type of field used to
                                store the moment and on the type of quadrature
                                node.
2015-05-23 Alberto Passalacqua: Added IOobject::groupname for improved naming
                                of files associated to moments.
2015-05-24 Alberto Passalacqua: Generalized moment update function to deal with
                                standard and extended nodes.
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

\*---------------------------------------------------------------------------*/

#include "moment.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template <class fieldType, class nodeType>
Foam::word
Foam::moment<fieldType, nodeType>::listToWord(const labelList& lst)
{
    word w;

    forAll(lst, dimi)
    {
        w += Foam::name(lst[dimi]);
    }

    return w;
}

template <class fieldType, class nodeType>
Foam::label
Foam::moment<fieldType, nodeType>::listToLabel(const labelList& lst)
{
    label l = 0;

    forAll(lst, dimi)
    {
        l += lst[dimi]*pow(10, lst.size() - dimi - 1);
    }

    return l;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class fieldType, class nodeType>
Foam::moment<fieldType, nodeType>::moment
(
    const word& distributionName,
    const labelList& cmptOrders,
    const fvMesh& mesh,
    const autoPtr<PtrList<nodeType>>& nodes
)
:
    fieldType
    (
        IOobject
        (
            momentName(listToWord(cmptOrders), distributionName),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    distributionName_(distributionName),
    nodes_(nodes),
    cmptOrders_(cmptOrders),
    name_(momentName(listToWord(cmptOrders_), distributionName_)),
    nDimensions_(cmptOrders_.size()),
    order_(sum(cmptOrders_))
{}

template <class fieldType, class nodeType>
Foam::moment<fieldType, nodeType>::moment
(
    const word& distributionName,
    const labelList& cmptOrders,
    const autoPtr<PtrList<nodeType>>& nodes,
    const fieldType& initMoment
)
:
    fieldType(initMoment),
    distributionName_(distributionName),
    nodes_(nodes),
    cmptOrders_(cmptOrders),
    name_(momentName(listToWord(cmptOrders_), distributionName_)),
    nDimensions_(cmptOrders_.size()),
    order_(sum(cmptOrders_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class fieldType, class nodeType>
Foam::moment<fieldType, nodeType>::~moment()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <class fieldType, class nodeType>
Foam::autoPtr<Foam::moment<fieldType, nodeType>>
Foam::moment<fieldType, nodeType>::clone() const
{
    NotImplemented;
    return autoPtr<moment<fieldType, nodeType>>(NULL);
}

template <class fieldType, class nodeType>
void Foam::moment<fieldType, nodeType>::update()
{
    // Resetting the moment to zero
    *this == dimensionedScalar("moment", (*this).dimensions(), 0);

    const PtrList<nodeType>& nodes = nodes_();

    bool extendedNode = nodes[0].extended();

    // If nodes do not have extended status, only use primary quadrature.
    if (!extendedNode)
    {
        forAll(nodes, pNodei)
        {
            const nodeType& node = nodes[pNodei];

            if (!node.extended())
            {
                fieldType m = node.primaryWeight();

                for (label cmpt = 0; cmpt < nDimensions_; cmpt++)
                {
                    const label cmptMomentOrder = cmptOrders()[cmpt];

                    tmp<fieldType> abscissaCmpt
                            = node.primaryAbscissa().component(cmpt);

                    tmp<fieldType> mPow = m*pow(abscissaCmpt, cmptMomentOrder);
                    m.dimensions().reset(mPow().dimensions());

                    m == mPow;
                }

                *this == *this + m;
            }
        }

        return;
    }

    // Extended quadrature case
    forAll(nodes, pNodei)
    {
        const nodeType& node = nodes[pNodei];

        const fieldType& pW = node.primaryWeight();

        for (label sNodei = 0; sNodei < node.nSecondaryNodes(); sNodei++)
        {
            fieldType m(pW*node.secondaryWeights()[sNodei]);

            for (label cmpt = 0; cmpt < nDimensions_; cmpt++)
            {
                const label cmptMomentOrder = cmptOrders()[cmpt];

                tmp<fieldType> abscissaCmpt
                        = node.secondaryAbscissae()[sNodei].component(cmpt);

                tmp<fieldType> mPow = m*pow(abscissaCmpt, cmptMomentOrder);

                m.dimensions().reset(mPow().dimensions());

                m == mPow;
            }

            *this == *this + m;
        }
    }
}


// ************************************************************************* //
