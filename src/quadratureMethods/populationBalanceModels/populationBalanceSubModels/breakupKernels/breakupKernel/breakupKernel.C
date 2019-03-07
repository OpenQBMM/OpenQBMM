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

\*---------------------------------------------------------------------------*/

#include "breakupKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(breakupKernel, 0);

    defineRunTimeSelectionTable(breakupKernel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernel::breakupKernel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    Cb_
    (
        dict.lookupOrDefault
        (
            "Cb",
            dimensionedScalar("one", inv(dimTime), 1.0)
        )
    ),
    daughterDistribution_
    (
        Foam::populationBalanceSubModels::daughterDistribution::New
        (
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::breakupKernel::~breakupKernel()
{}


Foam::scalar
Foam::populationBalanceSubModels::breakupKernel::breakupSource
(
    const label momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature
)
{
    scalar bSource = 0.0;

    const PtrList<volNode>& nodes = quadrature.nodes();

    if (!nodes[0].extended())
    {
        forAll(nodes, pNodeI)
        {
            const volNode& node = nodes[pNodeI];

            scalar bAbscissa = max(node.primaryAbscissae()[0][celli], 0.0);

            bSource += node.primaryWeight()[celli]
                    *Kb(bAbscissa, celli)
                    *(
                        daughterDistribution_->mD(momentOrder, bAbscissa)//Birth
                      - pow(bAbscissa, momentOrder)   //Death
                    );
        }

        return bSource;
    }

    forAll(nodes, pNodeI)
    {
        const volNode& node = nodes[pNodeI];

        forAll(node.secondaryWeights(), sNodei)
        {
            scalar bAbscissa
                = max(node.secondaryAbscissae()[0][sNodei][celli], 0.0);

            bSource += node.primaryWeight()[celli]
                *node.secondaryWeights()[sNodei][celli]
                *Kb(bAbscissa, celli)
                *(
                    daughterDistribution_->mD(momentOrder, bAbscissa)   //Birth
                  - pow(bAbscissa, momentOrder)                         //Death
                 );
        }
    }

    return bSource;
}

// ************************************************************************* //
