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
    Foam::velocityQuadratureNode

Description
    Implementation of the velocity quadrature node.

SourceFiles
    velocityQuadratureNode.C
    velocityQuadratureNodeI.H

\*---------------------------------------------------------------------------*/

#ifndef velocityQuadratureNode_H
#define velocityQuadratureNode_H

#include "quadratureNode.H"
#include "PtrList.H"
#include "dictionary.H"
#include "dictionaryEntry.H"
#include "mappedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class velocityQuadratureNode Declaration
\*---------------------------------------------------------------------------*/

template<class scalarType, class vectorType>
class velocityQuadratureNode
:
    public quadratureNode<scalarType, vectorType>
{
    // Private data

        //- Primary abscissa of the node
        vectorType velocityAbscissae_;


    // Private functions

        //- Lookup boundary conditions and return velocity abscissae field
        tmp<vectorType> createVelocityAbscissae
        (
            const scalarType& weight,
            const wordList& boundaryTypes = wordList()
        ) const;


public:

    // Constructors

        //- Construct from name, mesh and dimensions
        velocityQuadratureNode
        (
            const word& name,
            const word& distributionName,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const PtrList<dimensionSet>& abscissaDimensions,
            const wordList& boundaryTypes,
            const bool extended = false,
            const label nSecondaryNodes = 0
        );

        //- Construct from name, number of secondary nodes, mesh and dimensions
        velocityQuadratureNode
        (
            const word& name,
            const word& distributionName,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const PtrList<dimensionSet>& abscissaDimensions,
            const bool extended = false,
            const label nSecondaryNodes = 0
        );

        //- Return clone
        autoPtr<velocityQuadratureNode<scalarType, vectorType>> clone() const;

        //- Return a pointer to a new quadrature node created on freestore
        //  from Istream
        class iNew
        {
            const word distributionName_;
            const fvMesh& mesh_;
            const dimensionSet& weightDimensions_;
            const PtrList<dimensionSet>& abscissaeDimensions_;
            const wordList& boundaryTypes_;
            const bool extended_;
            const label nSecondaryNodes_;

        public:

            iNew
            (
                const word& distributionName,
                const fvMesh& mesh,
                const dimensionSet& weightDimensions,
                const PtrList<dimensionSet>& abscissaeDimensions,
                const wordList& boundaryTypes,
                const bool extended,
                const label nSecondaryNodes
            )
            :
                distributionName_(distributionName),
                mesh_(mesh),
                weightDimensions_(weightDimensions),
                abscissaeDimensions_(abscissaeDimensions),
                boundaryTypes_(boundaryTypes),
                extended_(extended),
                nSecondaryNodes_(nSecondaryNodes)
            {}

            autoPtr<velocityQuadratureNode<scalarType, vectorType>>
            operator()
            (
                Istream& is
            ) const
            {
                labelList ent(is);

                return autoPtr<velocityQuadratureNode<scalarType, vectorType>>
                (
                    new velocityQuadratureNode<scalarType, vectorType>
                    (
                        "node" + mappedList<scalar>::listToWord(ent),
                        distributionName_,
                        mesh_,
                        weightDimensions_,
                        abscissaeDimensions_,
                        boundaryTypes_,
                        extended_,
                        nSecondaryNodes_
                    )
                );
            }
        };


    //- Destructor
    virtual ~velocityQuadratureNode();


    // Member Functions

        // Access

            //- Const access to the velocity abscissa of the node
            virtual const vectorType& velocityAbscissae() const;

            //- Non-const access to the velcity abscissa of the node
            virtual vectorType& velocityAbscissae();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "velocityQuadratureNodeI.H"
#   include "velocityQuadratureNode.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
