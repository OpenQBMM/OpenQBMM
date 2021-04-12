/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2021 Alberto Passalacqua
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

#include "errorEstimator.H"
#include "staticFvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorEstimator::errorEstimator(const fvMesh& mesh)
:
    mesh_(mesh),
    needError_(!isA<staticFvMesh>(mesh)),
    dict_
    (
        needError_
      ? mesh.solutionDict().subDict("errorEstimate")
      : mesh.solutionDict()
    ),
    error_
    (
        IOobject
        (
            "error",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimless, 0.0)
    )
{
    if (needError_)
    {
        scalarFields_ =
            wordList
            (
                dict_.lookupOrDefault("scalarFields", wordList())
            );

        vectorFields_ =
            wordList
            (
                dict_.lookupOrDefault("vectorFields", wordList())
            );

        scalarScales_ =
            scalarField
            (
                dict_.lookupOrDefault
                ("scalarScaleFactors", scalarField(scalarFields_.size(), 1.0))
            );

        vectorScales_ =
            vectorField
            (
                dict_.lookupOrDefault
                (
                    "vectorScaleFactors",
                    vectorField(vectorFields_.size(), vector::one)
                )
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorEstimator::~errorEstimator()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::errorEstimator::estimateError()
{
    if (!needError_)
    {
        return;
    }

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    error_ = 0.0;

    forAll(scalarFields_, fieldi)
    {
        const volScalarField& f =
            mesh_.lookupObject<volScalarField>(scalarFields_[fieldi]);

        for (label facei = 0; facei < nInternalFaces; facei++)
        {
            label own = owner[facei];
            label nei = neighbour[facei];

            error_[own] =
                max
                (
                    error_[own],
                    mag(f[own] - f[nei])/scalarScales_[fieldi]
                );

            error_[nei] = max(error_[nei], error_[own]);
        }
    }
    forAll(vectorFields_, fieldi)
    {
        const volVectorField& f =
            mesh_.lookupObject<volVectorField>(vectorFields_[fieldi]);

        for (label facei = 0; facei < nInternalFaces; facei++)
        {
            label own = owner[facei];
            label nei = neighbour[facei];

            for (label cmpti = 0; cmpti < 3; cmpti++)
            {
                error_[own] =
                    max
                    (
                        error_[own],
                        mag
                        (
                            (f[own][cmpti] - f[nei][cmpti])
                           /vectorScales_[fieldi][cmpti]
                        )
                    );

                error_[nei] = max(error_[nei], error_[own]);
            }
        }
    }
}

// ************************************************************************* //
