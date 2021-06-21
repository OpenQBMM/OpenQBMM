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
    Copyright (C) 2019 Alberto Passalacqua
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
2017-03-26 Alberto Passalacqua: Added the capability to recompute the moment
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
        l += lst[dimi]*pow(scalar(10), lst.size() - dimi - 1);
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
    const autoPtr<mappedPtrList<nodeType>>& nodes
)
:
    fieldType
    (
        IOobject
        (
            momentName("moment", listToWord(cmptOrders), distributionName),
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
    name_(momentName("moment", listToWord(cmptOrders_), distributionName_)),
    nDimensions_(cmptOrders_.size()),
    order_(sum(cmptOrders_))
{}

template <class fieldType, class nodeType>
Foam::moment<fieldType, nodeType>::moment
(
    const word& distributionName,
    const labelList& cmptOrders,
    const autoPtr<mappedPtrList<nodeType>>& nodes,
    const fieldType& initMoment,
    const word momentSetName
)
:
    fieldType
    (
        momentName
        (
            "moment" + momentSetName,
            listToWord(cmptOrders),
            distributionName
        ),
        initMoment
    ),
    distributionName_(distributionName),
    nodes_(nodes),
    cmptOrders_(cmptOrders),
    name_
    (
        momentName
        (
            "moment" + momentSetName,
            listToWord(cmptOrders),
            distributionName
        )
    ),
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
    return nullptr;
}

template <class fieldType, class nodeType>
void Foam::moment<fieldType, nodeType>::update()
{
    // Resetting the moment to zero
    fieldType& moment(*this);
    moment == dimensionedScalar("moment", (*this).dimensions(), 0);

    const mappedPtrList<nodeType>& nodes = nodes_();
    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    // If nodes are not of extended type, only use primary quadrature.
    if (!nodes[0].extended())
    {
        forAll(nodes, pNodei)
        {
            const nodeType& node = nodes[pNodei];
            fieldType m = node.primaryWeight();

            for (label cmpt = 0; cmpt < scalarIndexes.size(); cmpt++)
            {
                label cmpti = scalarIndexes[cmpt];
                const label cmptMomentOrder = cmptOrders()[cmpti];

                tmp<fieldType> abscissaCmpt =
                    node.primaryAbscissae()[cmpt];

                tmp<fieldType> mPow = m*pow(abscissaCmpt, cmptMomentOrder);
                m.dimensions().reset(mPow().dimensions());

                m == mPow;
            }
            for (label cmpt = 0; cmpt < velocityIndexes.size(); cmpt++)
            {
                label cmpti = velocityIndexes[cmpt];
                const label cmptMomentOrder = cmptOrders()[cmpti];

                tmp<fieldType> abscissaCmpt
                        = node.velocityAbscissae().component(cmpt);

                tmp<fieldType> mPow = m*pow(abscissaCmpt, cmptMomentOrder);
                m.dimensions().reset(mPow().dimensions());

                m == mPow;
            }

            moment == moment + m;
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
            fieldType m(pW);

            for (label cmpt = 0; cmpt < scalarIndexes.size(); cmpt++)
            {
                label cmpti = scalarIndexes[cmpt];
                const label cmptMomentOrder = cmptOrders()[cmpti];

                tmp<fieldType> abscissaCmpt
                        = node.secondaryAbscissae()[cmpt][sNodei];

                tmp<fieldType> mPow =
                    m
                   *node.secondaryWeights()[cmpt][sNodei]
                   *pow(abscissaCmpt, cmptMomentOrder);

                m.dimensions().reset(mPow().dimensions());
                m == mPow;
            }
            for (label cmpt = 0; cmpt < velocityIndexes.size(); cmpt++)
            {
                label cmpti = velocityIndexes[cmpt];
                const label cmptMomentOrder = cmptOrders()[cmpti];

                tmp<fieldType> abscissaCmpt
                        = node.velocityAbscissae().component(cmpt);

                tmp<fieldType> mPow = m*pow(abscissaCmpt, cmptMomentOrder);
                m.dimensions().reset(mPow().dimensions());

                m == mPow;
            }
            moment == moment + m;
        }
    }
}

template <class fieldType, class nodeType>
void Foam::moment<fieldType, nodeType>::updateBoundaries()
{
    const mappedPtrList<nodeType>& nodes = nodes_();
    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    // If nodes are not of extended type, only use primary quadrature.
    if (!nodes[0].extended())
    {
        forAll(this->boundaryField(), patchi)
        {
            this->boundaryFieldRef()[patchi] = Zero;
            forAll(nodes, pNodei)
            {
                const nodeType& node = nodes[pNodei];
                scalarField m(node.primaryWeight().boundaryField()[patchi]);

                for (label cmpt = 0; cmpt < scalarIndexes.size(); cmpt++)
                {
                    label cmpti = scalarIndexes[cmpt];
                    const label cmptMomentOrder = cmptOrders()[cmpti];

                    const scalarField& abscissaCmpt =
                        node.primaryAbscissae()[cmpt].boundaryField()[patchi];

                    m *= pow(abscissaCmpt, cmptMomentOrder);

                }
                for (label cmpt = 0; cmpt < velocityIndexes.size(); cmpt++)
                {
                    label cmpti = velocityIndexes[cmpt];
                    const label cmptMomentOrder = cmptOrders()[cmpti];

                    tmp<scalarField> abscissaCmpt
                            = node.velocityAbscissae().boundaryField()[patchi].component(cmpt);

                    m *= pow(abscissaCmpt, cmptMomentOrder);
                }

                this->boundaryFieldRef()[patchi] += m;
            }
        }

        return;
    }

    // Extended quadrature case
    forAll(this->boundaryField(), patchi)
    {
        this->boundaryFieldRef()[patchi] = Zero;
        forAll(nodes, pNodei)
        {
            const nodeType& node = nodes[pNodei];
            const scalarField& pW = node.primaryWeight().boundaryField()[patchi];

            for (label sNodei = 0; sNodei < node.nSecondaryNodes(); sNodei++)
            {
                scalarField m(pW);

                for (label cmpt = 0; cmpt < scalarIndexes.size(); cmpt++)
                {
                    label cmpti = scalarIndexes[cmpt];
                    const label cmptMomentOrder = cmptOrders()[cmpti];

                    const scalarField& abscissaCmpt =
                        node.secondaryAbscissae()[cmpt][sNodei].boundaryField()[patchi];

                    m *=
                        node.secondaryWeights()[cmpt][sNodei].boundaryField()[patchi]
                       *pow(abscissaCmpt, cmptMomentOrder);
                }
                for (label cmpt = 0; cmpt < velocityIndexes.size(); cmpt++)
                {
                    label cmpti = velocityIndexes[cmpt];
                    const label cmptMomentOrder = cmptOrders()[cmpti];

                    tmp<scalarField> abscissaCmpt
                            = node.velocityAbscissae().boundaryField()[patchi].component(cmpt);

                    m *= pow(abscissaCmpt, cmptMomentOrder);

                }
                this->boundaryFieldRef()[patchi] += m;
            }
        }
    }
}

template <class fieldType, class nodeType>
void Foam::moment<fieldType, nodeType>::updateLocalMoment(label elemi)
{
    // Resetting the moment to zero
    scalar moment = 0;

    const mappedPtrList<nodeType>& nodes = nodes_();
    const labelList& scalarIndexes = nodes[0].scalarIndexes();
    const labelList& velocityIndexes = nodes[0].velocityIndexes();

    // If nodes are not of extended type, only use primary quadrature.
    if (!nodes[0].extended())
    {
        forAll(nodes, pNodei)
        {
            const nodeType& node = nodes[pNodei];
            scalar m = node.primaryWeight()[elemi];

            for (label cmpt = 0; cmpt < scalarIndexes.size(); cmpt++)
            {
                label cmpti = scalarIndexes[cmpt];
                const label cmptMomentOrder = cmptOrders()[cmpti];

                const scalar abscissaCmpt =
                    node.primaryAbscissae()[cmpt][elemi];

                m *= pow(abscissaCmpt, cmptMomentOrder);
            }
            for (label cmpt = 0; cmpt < velocityIndexes.size(); cmpt++)
            {
                label cmpti = velocityIndexes[cmpt];
                const label cmptMomentOrder = cmptOrders()[cmpti];

                const scalar abscissaCmpt =
                    component(node.velocityAbscissae()[elemi], cmpt);

                m *= pow(abscissaCmpt, cmptMomentOrder);
            }

            moment += m;
        }

        (*this)[elemi] = moment;

        return;
    }

    // Extended quadrature case
    forAll(nodes, pNodei)
    {
        const nodeType& node = nodes[pNodei];
        const scalar pW = node.primaryWeight()[elemi];

        for (label sNodei = 0; sNodei < node.nSecondaryNodes(); sNodei++)
        {
            scalar m = pW;

            for (label cmpt = 0; cmpt < nDimensions_; cmpt++)
            {
                label cmpti = scalarIndexes[cmpt];
                const label cmptMomentOrder = cmptOrders()[cmpti];

                const scalar abscissaCmpt =
                    node.secondaryAbscissae()[cmpt][sNodei][elemi];

                m *=
                    node.secondaryWeights()[cmpt][sNodei][elemi]
                   *pow(abscissaCmpt, cmptMomentOrder);
            }
            for (label cmpt = 0; cmpt < velocityIndexes.size(); cmpt++)
            {
                label cmpti = velocityIndexes[cmpt];
                const label cmptMomentOrder = cmptOrders()[cmpti];

                const scalar abscissaCmpt =
                    component(node.velocityAbscissae()[elemi], cmpt);

                m *= pow(abscissaCmpt, cmptMomentOrder);
            }
            moment += m;
        }
    }

    (*this)[elemi] = moment;
}

// ************************************************************************* //
