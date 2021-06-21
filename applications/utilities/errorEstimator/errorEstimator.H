/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::errorEstimator

Description
    Estimates error across faces with existing fields and scale factor for use
    with dynamicRefineFvMesh

SourceFiles
    errorEstimator.C

\*---------------------------------------------------------------------------*/

#ifndef errorEstimator_H
#define errorEstimator_H

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class errorEstimator Declaration
\*---------------------------------------------------------------------------*/

class errorEstimator
{
    // Private data

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Is the mesh moving/refining
        bool needError_;

        //- Reference to error dictionary
        const dictionary& dict_;

        //- Stored error field
        volScalarField error_;

        //- List of scalar field names
        wordList scalarFields_;

        //- List of vector field names
        wordList vectorFields_;

        //- List of scalar scale factors
        scalarField scalarScales_;

        //- List of vector scale factors
        vectorField vectorScales_;


public:

    // Constructors

        //- Construct from a given vector
        errorEstimator(const fvMesh& mesh);


    //- Destructor
    ~errorEstimator();


    // Member Functions

        //- Solve the errorEstimator linear system with the given source vector
        void estimateError();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
