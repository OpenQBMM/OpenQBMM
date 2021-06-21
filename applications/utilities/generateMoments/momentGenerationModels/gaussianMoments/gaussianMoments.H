/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 Alberto Passalacqua
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
    Foam::momentGenerationSubModels::gaussian

Description
    Construct moments using mean velocity and covariance tensor

SourceFiles
    gaussianMoments.C

\*---------------------------------------------------------------------------*/

#ifndef gaussian_H
#define gaussian_H

#include "momentGenerationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{

/*---------------------------------------------------------------------------*\
                        Class gaussian Declaration
\*---------------------------------------------------------------------------*/

class gaussian
:
    public momentGenerationModel
{
protected:
    // Protected data

        //- Zero order moment
        volScalarField m0_;

        //- Mean velocity field
        volVectorField U_;

        //- Isotropic variance
        volScalarField Theta_;

        //- Full covariance tensor
        volSymmTensorField Sigma_;

        //- Flag if isotropic
        bool isotropic_;

        //- Number of velocity directions
        label nVelocityDimensions_;

        //- Set nodes give alpha, U, Thea, and Sigma
        void setNodes
        (
            const dictionary& dict,
            const scalarField& alpha,
            const vectorField& U,
            const scalarField& Theta,
            const symmTensorField& Sigma
        );

public:

        //- Runtime type information
        TypeName("gaussian");


    // Constructors

        //- Construct from components
        gaussian
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const labelListList& momentOrders,
            const label nNodes
        );


    //- Destructor
    virtual ~gaussian();


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
