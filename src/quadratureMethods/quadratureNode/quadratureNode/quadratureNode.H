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
#include "mappedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class quadratureNode Declaration
\*---------------------------------------------------------------------------*/
template<class scalarType, class vectorType>
class quadratureNode
{
protected:

        //- Typedefs

        typedef scalarType weightType;
        typedef PtrList<PtrList<weightType>> secondaryWeightsType;

        typedef scalarType abscissaType;
        typedef PtrList<abscissaType> abscissaeType;
        typedef PtrList<abscissaeType> secondaryAbscissaeType;

        typedef vectorType velocityAbscissaeType;

        typedef scalarType sigmaType;
        typedef PtrList<sigmaType> sigmasType;


    // Protected data

        //- Name of the quadrature node
        word name_;

        //- Primary weight of the node
        weightType weight_;

        //- Primary abscissae of the node
        abscissaeType abscissae_;

        //- Indicies of abscissa components
        labelList scalarIndexes_;

        //- Indices of velocity components
        labelList velocityIndexes_;

        //- Index of size component
        label sizeIndex_;

        //- Bool to indicate if abscissas are based on length
        bool lengthBased_;

        //- Bool to indicate if abscissas are based on mass
        bool massBased_;

        //- Pointer to density field
        const volScalarField* rhoPtr_;

        //- Bool to indicate if weights are volume fractions
        bool useVolumeFraction_;

        //- List of secondary weights of the node
        secondaryWeightsType secondaryWeights_;

        //- List of secondary abscissae of the node
        secondaryAbscissaeType secondaryAbscissae_;

        //- Parameter of the kernel density function used in EQMOM
        sigmasType sigmas_;

        //- Number of secondary nodes
        label nSecondaryNodes_;

        //- Boolean variable to identify extended nodes
        bool extended_;



public:

    // Constructors

        // Constructors

        //- Construct from name, mesh and dimensions
        quadratureNode
        (
            const word& name,
            const word& distributionName,
            const fvMesh& mesh,
            const dimensionSet& weightDimensions,
            const PtrList<dimensionSet>& abscissaeDimensions,
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
            const PtrList<dimensionSet>& abscissaeDimensions,
            const bool extended = false,
            const label nSecondaryNodes = 0
        );

        //- Return clone
        autoPtr<quadratureNode<scalarType, vectorType>> clone() const;

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

            autoPtr<quadratureNode<scalarType, vectorType>>
            operator()
            (
                Istream& is
            ) const
            {
                labelList ent(is);
                
                return autoPtr<quadratureNode<scalarType, vectorType>>
                (
                    new quadratureNode<scalarType, vectorType>
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
    virtual ~quadratureNode();


    // Member Functions

        // Access

            //- Return the name of the node
            const word& name() const;

            //- Return abscissa indexes
            const labelList& scalarIndexes() const;

            //- Return velocity indexes
            const labelList& velocityIndexes() const;

            //- Return size index
            label sizeIndex() const;

            //- Is the node length based
            bool lengthBased() const;

            //- Is the weight in terms of volume fraction
            bool useVolumeFraction() const;

            //- Is the node extended
            bool extended() const;

            //- Returns the number of secondary nodes
            label nSecondaryNodes() const;

            //- Const access to the weight of the node
            const weightType& primaryWeight() const;

            //- Non-const access to the weight of the node
            weightType& primaryWeight();

            //- Const access to the scalar abscissa of the node
            const abscissaeType& primaryAbscissae() const;

            //- Non-const access to the scalar abscissa of the node
            abscissaeType& primaryAbscissae();

            //- Const access to the secondary weights of the node
            const secondaryWeightsType& secondaryWeights() const;

            //- Non-const access to the secondary weights of the node
            secondaryWeightsType& secondaryWeights();

            //- Const access to the list of secondary abscissae of the node
            const secondaryAbscissaeType& secondaryAbscissae() const;

            //- Non-const access to the list of secondary abscissae of the node
            secondaryAbscissaeType& secondaryAbscissae();

            //- Const access to sigma
            const sigmasType& sigmas() const;

            //- Non-const access to sigma
            sigmasType& sigmas();

            //- Const access to the velocity abscissa of the node
            virtual const velocityAbscissaeType& velocityAbscissae() const;

            //- Non-const access to the velcity abscissa of the node
            virtual velocityAbscissaeType& velocityAbscissae();

            //- Return the diameter given an abscissae
            tmp<volScalarField> d(const volScalarField& abscissa) const;

            //- Return the diameter given an abscissae and celli
            scalar d(const label celli, const scalar& abscissa) const;

            //- Return the number density given an abscissae
            tmp<volScalarField> n
            (
                const volScalarField& weight,
                const volScalarField& abscissa
            ) const;

            //- Return the number density given an abscissae and celli
            scalar n
            (
                const label celli, 
                const scalar& weight, 
                const scalar& abscissa
            ) const;
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
