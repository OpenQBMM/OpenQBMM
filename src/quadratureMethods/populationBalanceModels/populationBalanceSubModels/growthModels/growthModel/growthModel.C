/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Matteo Icardi
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

\*---------------------------------------------------------------------------*/

#include "growthModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(growthModel, 0);

    defineRunTimeSelectionTable(growthModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModel::growthModel
(
    const dictionary& dict
)
:
    dict_(dict),
    Cg_
    (
        dict.lookupOrDefault
        (
            "Cg",
            dimensionedScalar("one", inv(dimTime), 1.0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::growthModel::~growthModel()
{}


Foam::scalar
Foam::populationBalanceSubModels::growthModel::phaseSpaceConvection
(
    const label momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    scalar gSource = 0.0;

    if (momentOrder < 1)
    {
        return gSource;
    }

    const PtrList<volNode>& nodes = quadrature.nodes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volNode& node = nodes[pNodeI];

            scalar bAbscissa = max(node.primaryAbscissae()[0][celli], 0.0);

            gSource += node.primaryWeight()[celli]
                    *Kg(node.primaryAbscissae()[0][celli])
                    *momentOrder*pow(bAbscissa, momentOrder - 1);
        }

        return gSource;
    }

    forAll(nodes, pNodeI)
    {
        const volNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            scalar bAbscissa
                = max(node.secondaryAbscissae()[0][sNodei][celli], 0.0);

            gSource += node.primaryWeight()[celli]
                *node.secondaryWeights()[sNodei][celli]
                *Kg(bAbscissa)
                *momentOrder*pow
                    (
                        bAbscissa,
                        momentOrder - 1
                    );
        }
    }

    return gSource;
}

// ************************************************************************* //
