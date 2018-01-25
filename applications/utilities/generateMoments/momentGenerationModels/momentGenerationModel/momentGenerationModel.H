/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 Alberto Passalacqua
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
    Foam::momentGenerationModel

Description
    Base class for runtime selection for moment generation. One object is used
    for all node. Moments are constructed outside of class.

SourceFiles
    momentGenetationModel.C
    momentGenerationModelnew.C

\*---------------------------------------------------------------------------*/

#ifndef momentGeneration_H
#define momentGeneration_H

#include "dictionary.H"
#include "mappedLists.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class momentGenerationModel
{

    // Private member functions

        //- Disallow default bitwise copy construct
        momentGenerationModel(const momentGenerationModel&);

        //- Disallow default bitwise assignment
        void operator=(const momentGenerationModel&);

        //- Create map for velocity nodes
        labelListList createNodeIndexes() const;


protected:

    // Protected data

        //- Constant refrence to dictionary
        const dictionary& dict_;

        //- Number of dimensions in distribution
        const label nDims_;

        //- Number of nodes
        const label nNodes_;

        //- Number of moments
        const label nMoments_;

        //- Node indexes
        labelListList nodeIndexes_;

        //- Moment orders
        labelListList momentOrders_;

        //- List of weights
        List<scalar> weights_;

        //- List of abscissae
        List<scalarList> abscissae_;

        //- List of moments
        mappedList<dimensionedScalar> moments_;


    // Protected member functions

        //- Update moments
        void updateMoments();

        //- Zero weights, abscissae and moments
        void reset();


public:

    //- Runtime information
    TypeName("momentGenerationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        momentGenerationModel,
        dictionary,
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const label nNodes
        ),
        (dict, momentOrders, nNodes)
    );


    // Constructor

        //- Construct from dictionary and number of nodes
        momentGenerationModel
        (
            const dictionary& dict,
            const labelListList& momentsOrders,
            const label nNodes
        );


    // Selector

        static autoPtr<momentGenerationModel> New
        (
            const dictionary& dict,
            const labelListList& momentsOrders,
            const label nNodes
        );


    //- Destructor
    virtual ~momentGenerationModel();


    // Public Functions

        // Access

            //- Return constant access to moments
            const mappedList<dimensionedScalar>& moments() const
            {
                return moments_;
            }

        // Member Functions

            //- Update weights and abscissae given a dictionary
            virtual void updateQuadrature(const dictionary&) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
