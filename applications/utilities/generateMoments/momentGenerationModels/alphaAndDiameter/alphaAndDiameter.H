/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 Alberto Passalacqua
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
    Foam::momentGenerationSubModels::alphaAndDiameter

Description
    Moment generation done using volume fraction, density, and bubble diameter.

SourceFiles
    alphaAndDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef alphaAndDiameter_H
#define alphaAndDiameter_H

#include "momentGenerationModel.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{

/*---------------------------------------------------------------------------*\
                    Class alphaAndDiameter Declaration
\*---------------------------------------------------------------------------*/

class alphaAndDiameter
:
    public momentGenerationModel
{
protected:

    //- Read alpha field
    volScalarField alpha_;

    //- Are nodes scaled by the total volume fraction
    Switch scale_;

    //- Copy of rho field
    volScalarField rho_;

    //- List of diameters
    List<scalarField> diameters_;

    //- List of volume fractions
    List<scalarField> alphas_;

    //- Total node volume fraction
    scalarField sumAlpha_;

    //- Is the moment set mass based
    Switch massBased_;


public:

        //- Runtime type information
        TypeName("alphaAndDiameter");


    // Constructors

        //- Construct from components
        alphaAndDiameter
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const labelListList& momentsOrders,
            const label nNodes
        );


    //- Destructor
    virtual ~alphaAndDiameter();


    // Public member functions

        //- Update weights and abscissae given a dictionary
        virtual void updateMoments
        (
            const dictionary& dict,
            const label patchi = -1
        );

        //- Update weights and abscissae given a dictionary and cell set
        virtual void updateMoments
        (
            const dictionary& dict,
            const labelList& cells
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace momentGenerationSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
