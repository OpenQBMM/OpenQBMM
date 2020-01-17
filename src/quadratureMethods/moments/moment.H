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

Class
    Foam::moment

Description
    Stores a moment of a distribution function with arbitrary dimensionality.

SourceFiles
    momentI.H
    moment.C

\*---------------------------------------------------------------------------*/

#ifndef moment_H
#define moment_H

#include "mappedPtrList.H"
#include "quadratureNodes.H"
#include "volFields.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class moment Declaration
\*---------------------------------------------------------------------------*/

template <class fieldType, class nodeType>
class moment
:
    public fieldType
{
    //- Name of the distribution associated to the moment
    const word distributionName_;

    //- Reference to the quadrature nodes
    const autoPtr<mappedPtrList<nodeType>>& nodes_;

    //- Order of each dimension of the moment
    const labelList cmptOrders_;

    //- Name of the moment
    const word name_;

    //- Number of dimensions of the distribution function
    const label nDimensions_;

    //- Overall order of the moment
    const label order_;


public:

    //- Convert a list of labels to a word. {1, 2, 3} -> 123
    static word listToWord(const labelList& lst);

    //- Convert a list of labels to a single label. {0, 2, 3} -> 23
    static label listToLabel(const labelList& lst);


    // Constructors

        //- Construct from a list of orders of each dimension of the moment
        //  and the quadrature nodes
        moment
        (
            const word& distributionName,
            const labelList& cmptOrders,
            const fvMesh& mesh,
            const autoPtr<mappedPtrList<nodeType>>& nodes
        );

        //- Construct from name, orders, and quadrature nodes
        moment
        (
            const word& distributionName,
            const labelList& cmptOrders,
            const autoPtr<mappedPtrList<nodeType>>& nodes,
            const fieldType& initMoment,
            const word momentSetName = word::null
        );

        //- Return clone
        autoPtr<moment> clone() const;


        //- Return a pointer to a new moment created on freestore
        //  from Istream
        class iNew
        {
            const word& distributionName_;
            const fvMesh& mesh_;
            const autoPtr<mappedPtrList<nodeType>>& nodes_;

        public:

            iNew
            (
                const word& distributionName,
                const fvMesh& mesh,
                const autoPtr<mappedPtrList<nodeType>>& nodes
            )
            :
                distributionName_(distributionName),
                mesh_(mesh),
                nodes_(nodes)
            {}

            autoPtr<moment<fieldType, nodeType>> operator()(Istream& is) const
            {
                labelList ent(is);
                return autoPtr<moment<fieldType, nodeType>>
                (
                    new moment<fieldType, nodeType>
                    (
                        distributionName_,
                        ent,
                        mesh_,
                        nodes_
                    )
                );
            }
        };


    //- Destructor
    ~moment();


    // Member Functions

            //- Creates the name of the moment
            word momentName
            (
                const word& momentSetName,
                const word& order,
                const word& distributionName
            )
            {
                word tmpName = IOobject::groupName(momentSetName, order);
                return IOobject::groupName(tmpName, distributionName);
            }

        // Access

            //- Return the name of the moment
            inline const word& name() const;

            //- Return the name of the moment
            inline const word& keyword() const;

            //- Return the dimensionality of the distribution
            inline label nDimensions() const;

            //- Return the overall order of the moment
            inline label order() const;

            //- Return the order of each dimension
            inline const labelList& cmptOrders() const;


        // Edit

            //- Update the moment from the quadrature nodes
            void update();

            //- Update the moment from the quadrature nodes
            //  on boundaries
            void updateBoundaries();

            //- Update the moment from the quadrature nodes in a single element
            void updateLocalMoment(label elemi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "momentI.H"
#   include "moment.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
