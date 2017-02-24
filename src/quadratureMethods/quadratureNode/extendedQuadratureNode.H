/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
    Foam::extendedQuadratureNode

Description
    Extended version of the quadrature node, which stores the list of secondary
    weights and abscissae associated to the primary weight and abscissa, and
    the sigma parameter.

SourceFiles
    extendedQuadratureNode.C
    extendedQuadratureNodeI.H

\*---------------------------------------------------------------------------*/

#ifndef extendedQuadratureNode_H
#define extendedQuadratureNode_H

#include "fvMesh.H"
#include "PtrList.H"
#include "dictionary.H"
#include "dictionaryEntry.H"
#include "quadratureNode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class extendedQuadratureNode Declaration
\*---------------------------------------------------------------------------*/
template <class weightType, class abscissaType, class sigmaType>
class extendedQuadratureNode
:
    public quadratureNode<weightType, abscissaType>
{
    // Private data

        //- Name of the quadrature node
        word name_;

        //- Number of secondary nodes
        label nSecondaryNodes_;

        //- List of secondary weights of the node
        PtrList<weightType> secondaryWeights_;

        //- List of secondary abscissae of the node
        PtrList<abscissaType> secondaryAbscissae_;

        //- Parameter of the kernel density function used in EQMOM
        sigmaType sigma_;

        //- Boolean variable to identify extended nodes
        bool extended_;

public:

    // Constructors

        //- Construct from name, dictionary, mesh and dimensions
        extendedQuadratureNode
        (
            const word& name,
            const word& distributionName,
            const dictionary& nodeDict,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const dimensionSet& abscissaDimensions,
            const wordList& boundaryTypes
        );

        //- Construct from name, number of secondary nodes, mesh and dimensions
        extendedQuadratureNode
        (
            const word& name,
            const word& distributionName,
            const label nSecondaryNodes,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const dimensionSet& abscissaDimensions
        );

        //- Return clone
        autoPtr<extendedQuadratureNode<weightType, abscissaType, sigmaType>>
            clone() const;

        //- Return a pointer to a new quadrature node created on freestore
        //  from Istream
        class iNew
        {
            const word distributionName_;
            const fvMesh& mesh_;
            const dimensionSet& weightDimensions_;
            const dimensionSet& abscissaDimensions_;
            const wordList& boundaryTypes_;

        public:

            iNew
            (
                const word& distributionName,
                const fvMesh& mesh,
                const dimensionSet& weightDimensions,
                const dimensionSet& abscissaDimensions,
                const wordList& boundaryTypes
            )
            :
                distributionName_(distributionName),
                mesh_(mesh),
                weightDimensions_(weightDimensions),
                abscissaDimensions_(abscissaDimensions),
                boundaryTypes_(boundaryTypes)
            {}

            autoPtr
            <
                extendedQuadratureNode<weightType, abscissaType, sigmaType>
            >
            operator()
            (
                Istream& is
            ) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr
                <
                    extendedQuadratureNode<weightType, abscissaType, sigmaType>
                >
                (
                    new extendedQuadratureNode
                    <
                        weightType, abscissaType, sigmaType
                    >
                    (
                        ent.keyword(),
                        distributionName_,
                        ent,
                        mesh_,
                        weightDimensions_,
                        abscissaDimensions_,
                        boundaryTypes_
                    )
                );
            }
        };


    //- Destructor
    virtual ~extendedQuadratureNode();


    // Member Functions

        // Access

            //- Return the name of the node
            inline const word& name() const;

            inline bool extended() const;

            //- Returns the number of secondary nodes
            inline label nSecondaryNodes() const;

            //- Const access to the secondary weights of the node
            inline const PtrList<weightType>& secondaryWeights() const;

            //- Non-const access to the secondary weights of the node
            inline PtrList<weightType>& secondaryWeights();

            //- Const access to the list of secondary abscissae of the node
            inline const PtrList<abscissaType>& secondaryAbscissae() const;

            //- Non-const access to the list of secondary abscissae of the node
            inline PtrList<abscissaType>& secondaryAbscissae();

            //- Const access to sigma
            inline const sigmaType& sigma() const;

            //- Non-const access to sigma
            inline sigmaType& sigma();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "extendedQuadratureNodeI.H"
#   include "extendedQuadratureNode.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
