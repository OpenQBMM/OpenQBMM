/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2015-02-19 Alberto Passalacqua: Templated class on type of weight and abscissa.
2015-03-08 Alberto Passalacqua: Generalised implementation to include secondary
                                quadrature and the sigma parameter of kernel 
                                density functions used in the extended 
                                quadrature method of moments.
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
    Stores the primary weight and abscissa, the list of secondary weights
    and abscissae associated to the primary weight and abscissa, and the sigma 
    parameter of an extended quadrature node.

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
        word name_;

        //- Number of secondary nodes
        label nSecondaryNodes_;

        //- Primary weight of the node
        weightType primaryWeight_;

        //- Primary abscissa of the node
        abscissaType primaryAbscissa_;

        //- List of secondary weights of the node
        PtrList<weightType> secondaryWeights_;

        //- List of secondary abscissae of the node
        PtrList<abscissaType> secondaryAbscissae_;

        //- Parameter of the kernel density function used in EQMOM
        //  It is assumed that the sigma parameter has the same type of 
        //  the quadrature weight. 
        weightType sigma_;

        //- Boolean variable to identify extended nodes
        bool extended_;

public:

    // Constructors

        //- Construct from name, dictionary, mesh and dimensions
        quadratureNode
        (
            const word name,
            const dictionary nodeDict,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const dimensionSet& abscissaDimensions
        );
        
        quadratureNode
        (
            const word name,
            const label nSecondaryNodes,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const dimensionSet& abscissaDimensions
        );

        //- Return clone
        autoPtr<quadratureNode<weightType, abscissaType> > clone() const;
        
        //- Return a pointer to a new quadrature node created on freestore
        //  from Istream
        class iNew
        {
            const fvMesh& mesh_;
            const dimensionSet& weightDimensions_;
            const dimensionSet& abscissaDimensions_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const dimensionSet& weightDimensions,
                const dimensionSet& abscissaDimensions
            )
            :
                mesh_(mesh),
                weightDimensions_(weightDimensions),
                abscissaDimensions_(abscissaDimensions)
            {}

            autoPtr<quadratureNode<weightType, abscissaType> > operator()
            (
                Istream& is
            ) const
            {
                dictionaryEntry ent(dictionary::null, is);
                return autoPtr<quadratureNode<weightType, abscissaType> >
                (
                    new quadratureNode<weightType, abscissaType>
                    (
                        ent.keyword(), 
                        ent, 
                        mesh_,
                        weightDimensions_,
                        abscissaDimensions_
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

            //- Const access to the primary weight of the node
            inline const weightType& primaryWeight() const;

            //- Non-const access to the primary weight of the node
            inline weightType& primaryWeight();

            //- Const access to the primary abscissa of the node
            inline const abscissaType& primaryAbscissa() const;

            //- Non-const access to the primary abscissa of the node
            inline abscissaType& primaryAbscissa();

            //- Const access to the secondary weight of the node
            inline const PtrList<weightType>& secondaryWeights() const;

            //- Non-const access to the secondary weight of the node
            inline PtrList<weightType>& secondaryWeights();

            //- Const access to the list of secondary abscissae of the node
            inline const PtrList<abscissaType>& secondaryAbscissae() const;

            //- Non-const access to the list of secondary abscissae of the node
            inline PtrList<abscissaType>& secondaryAbscissae();

            //- Const access to sigma
            inline const weightType& sigma() const;

            //- Non-const access to sigma
            inline weightType& sigma();
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
