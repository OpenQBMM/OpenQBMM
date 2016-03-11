/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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
    Foam::mixingModel

Description
    Solve a turbulent mixing model.

SourceFiles
    mixingModel.C

\*---------------------------------------------------------------------------*/

#ifndef mixingModel_H
#define mixingModel_H

#include "univariatePDFTransportModel.H"
#include "mixingKernel.H"
#include "diffusionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace mixingModels
{

/*---------------------------------------------------------------------------*\
                        Class mixingModel Declaration
\*---------------------------------------------------------------------------*/

class mixingModel
:
    public univariatePDFTransportModel
{
    // Private data

        const word name_;

        //- Aggregation kernel model
        autoPtr
        <
            Foam::mixingSubModels::mixingKernel
        > mixingKernel_;

        //- Diffusion model
        autoPtr
        <
            Foam::mixingSubModels::diffusionModel
        > diffusionModel_;

    // Private member functions

        //- Diffusion term
        virtual tmp<fvScalarMatrix> momentDiffusion
        (
            const volUnivariateMoment& moment
        );

        //- Convect moments in phase space
        virtual tmp<volScalarField> phaseSpaceConvection
        (
            const volUnivariateMoment& moment
        );

        //- Calculate source terms
        virtual tmp<fvScalarMatrix> momentSource
        (
            const volUnivariateMoment& moment
        );


public:
    // Constructors

        //- Construct from components
        mixingModel
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~mixingModel();

    // Member Functions

        //- Solve population balance equation
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceModels
} // End namespace univariatePDFTransportModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
