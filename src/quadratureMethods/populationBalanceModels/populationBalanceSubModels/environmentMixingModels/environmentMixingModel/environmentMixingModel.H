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
    Copyright (C) 2019 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::environmentMixingModel

Description
    Abstract class for mixing model kernels.

SourceFiles
    environmentMixingModel.C

\*---------------------------------------------------------------------------*/

#ifndef environmentMixingModel_H
#define environmentMixingModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "momentFieldSets.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{

/*---------------------------------------------------------------------------*\
                    Class environmentMixingModel Declaration
\*---------------------------------------------------------------------------*/

class environmentMixingModel
{
protected:

    // Protected data

        //- Reference to dictionary
        const dictionary& dict_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Coefficient of aggregation kernel
        const dimensionedScalar Cphi_;


public:

    //- Runtime type information
    TypeName("environmentMixingModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        environmentMixingModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh
        ),
        (dict, mesh)
    );


    // Constructors

        //- Construct from components
        environmentMixingModel
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        environmentMixingModel(const environmentMixingModel&) = delete;


    // Selectors

        static autoPtr<environmentMixingModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~environmentMixingModel();


    // Member Functions

        //- Mixing kernel
        virtual tmp<fvScalarMatrix> K
        (
            const volScalarField& meanMoment,
            const volScalarField& meanMomentVariance,
            const volScalarField& meanMixtureFraction
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const environmentMixingModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
