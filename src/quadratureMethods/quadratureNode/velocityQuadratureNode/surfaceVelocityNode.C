/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/

#include "surfaceVelocityNode.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
tmp<surfaceVectorField>
velocityQuadratureNode<surfaceScalarField, surfaceVectorField>::
createVelocityAbscissae
(
    const surfaceScalarField& weight,
    const wordList& boundaryTypes
) const
{
    const fvMesh& mesh = weight.mesh();

    if (boundaryTypes.size() == 0)
    {
        return tmp<surfaceVectorField>
        (
            new surfaceVectorField
            (
                IOobject
                (
                    IOobject::groupName("velocityAbscissae", this->name_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector
                (
                    "zeroVelocityAbscissa",
                    dimVelocity,
                    Zero
                )
            )
        );
    }

    return tmp<surfaceVectorField>
    (
        new surfaceVectorField
        (
            IOobject
            (
                IOobject::groupName("velocityAbscissae", this->name_),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "zeroVelocityAbscissa",
                dimVelocity,
                Zero
            ),
            boundaryTypes
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
