/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2015-02-19 Alberto Passalacqua: Templated class on type of weight and abscissa
                                to store a basic quadrature node made of weight
                                and abscissa.
2015-11-02 Alberto Passalacqua: Generalized initialization of fields based on
                                value_type of the field.
2015-11-05 Alberto Passalacqua: Added interface functions for secondary
                                quadrature (not implemented).
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
    Stores a quadrature weight and a quadrature abscissa, whose type is
    templated.

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
template <class weightType, class abscissaType>
class quadratureNode
{
    // Private data

        //- Name of node
        const word name_;

        //- Node dictionary
        dictionary nodeDict_;

        //- Primary weight of the node
        weightType weight_;

        //- Primary abscissa of the node
        abscissaType abscissa_;

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
            const wordList& boundaryTypes
        );

        //- Construct from name, mesh and dimensions
        quadratureNode
        (
            const word& name,
            const word& distributionName,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const dimensionSet& abscissaDimensions
        );

        //- Return clone
        autoPtr<quadratureNode<weightType, abscissaType>> clone() const;

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

            autoPtr<quadratureNode<weightType, abscissaType>> operator()
            (
                Istream& is
            ) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<quadratureNode<weightType, abscissaType>>
                (
                    new quadratureNode<weightType, abscissaType>
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

            //- Const access to the secondary weight of the node
            inline const PtrList<weightType>& secondaryWeights() const;

            //- Non-const access to the secondary weight of the node
            inline PtrList<weightType>& secondaryWeights();

            //- Const access to the list of secondary abscissae of the node
            inline const PtrList<abscissaType>& secondaryAbscissae() const;

            //- Non-const access to the list of secondary abscissae of the node
            inline PtrList<abscissaType>& secondaryAbscissae();
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
