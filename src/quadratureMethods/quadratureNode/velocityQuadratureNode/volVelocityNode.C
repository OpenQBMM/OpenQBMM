/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

InClass

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#include "volVelocityNode.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
tmp<volVectorField>
velocityQuadratureNode<volScalarField, volVectorField>::
createVelocityAbscissae
(
    const volScalarField& weight,
    const wordList& boundaryTypes
) const
{
    const fvMesh& mesh = weight.mesh();
    word UName = IOobject::groupName("U", weight.group());

    if (mesh.foundObject<volVectorField>(UName) && boundaryTypes.size() > 0)
    {
        const volVectorField& Umean =
            mesh.lookupObject<volVectorField>(UName);

        return tmp<volVectorField>
        (
            new volVectorField
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
                Umean.boundaryField().types()
            )
        );
    }

    return tmp<volVectorField>
    (
        new volVectorField
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
