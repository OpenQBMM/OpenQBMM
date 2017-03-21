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
    Foam::quadratureNode

Description
    Extended version of the quadrature node, which stores the list of secondary
    weights and abscissae associated to the primary weight and abscissa, and
    the sigma parameter.

SourceFiles
    quadratureNode.C
    quadratureNodeI.H

\*---------------------------------------------------------------------------*/

#ifndef quadratureNode_H
#define quadratureNode_H

#include "fvMesh.H"
#include "PtrList.H"
#include "dictionary.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class quadratureNode Declaration
\*---------------------------------------------------------------------------*/
template <class weightType, class abscissaType, class sigmaType>
class quadratureNode
{
    // Private data

        //- Name of the quadrature node
        word name_;

        //- Node dictionary
        dictionary nodeDict_;

        //- Primary weight of the node
        weightType weight_;

        //- Primary abscissa of the node
        abscissaType abscissa_;

        //- List of secondary weights of the node
        PtrList<weightType> secondaryWeights_;

        //- List of secondary abscissae of the node
        PtrList<abscissaType> secondaryAbscissae_;

        //- Parameter of the kernel density function used in EQMOM
        autoPtr<sigmaType> sigma_;

        //- Number of secondary nodes
        label nSecondaryNodes_;

        //- Boolean variable to identify extended nodes
        bool extended_;

public:

    // Constructors

        //- Construct from name, dictionary, mesh and dimensions
        quadratureNode
        (
            const word& name,
            const word& distributionName,
            const dictionary& nodeDict,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const dimensionSet& abscissaDimensions,
            const wordList& boundaryTypes,
            const bool extended = false,
            const label nSecondaryNodes = 0
        );

        //- Construct from name, number of secondary nodes, mesh and dimensions
        quadratureNode
        (
            const word& name,
            const word& distributionName,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const dimensionSet& abscissaDimensions,
            const bool extended = false,
            const label nSecondaryNodes = 0
        );

        //- Return clone
        autoPtr<quadratureNode<weightType, abscissaType, sigmaType>>
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
            const bool extended_;
            const label nSecondaryNodes_;

        public:

            iNew
            (
                const word& distributionName,
                const fvMesh& mesh,
                const dimensionSet& weightDimensions,
                const dimensionSet& abscissaDimensions,
                const wordList& boundaryTypes,
                const bool extended,
                const label nSecondaryNodes
            )
            :
                distributionName_(distributionName),
                mesh_(mesh),
                weightDimensions_(weightDimensions),
                abscissaDimensions_(abscissaDimensions),
                boundaryTypes_(boundaryTypes),
                extended_(extended),
                nSecondaryNodes_(nSecondaryNodes)
            {}

            autoPtr
            <
                quadratureNode<weightType, abscissaType, sigmaType>
            >
            operator()
            (
                Istream& is
            ) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr
                <
                    quadratureNode<weightType, abscissaType, sigmaType>
                >
                (
                    new quadratureNode
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
                        boundaryTypes_,
                        extended_,
                        nSecondaryNodes_
                    )
                );
            }
        };


    //- Destructor
    virtual ~quadratureNode();


    // Member Functions

        // Access

            //- Return the name of the node
            inline const word& name() const;

            inline bool extended() const;

            //- Returns the number of secondary nodes
            inline label nSecondaryNodes() const;

            //- Const access to the weight of the node
            inline const weightType& primaryWeight() const;

            //- Non-const access to the weight of the node
            inline weightType& primaryWeight();

            //- Const access to the abscissa of the node
            inline const abscissaType& primaryAbscissa() const;

            //- Non-const access to the abscissa of the node
            inline abscissaType& primaryAbscissa();

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
#   include "quadratureNodeI.H"
#   include "quadratureNode.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
