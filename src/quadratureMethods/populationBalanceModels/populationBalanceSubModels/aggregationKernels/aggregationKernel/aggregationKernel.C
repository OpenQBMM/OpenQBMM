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

#include "aggregationKernel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
    defineTypeNameAndDebug(aggregationKernel, 0);

    defineRunTimeSelectionTable(aggregationKernel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernel::aggregationKernel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh),
    Ca_
    (
        dict.lookupOrDefault
        (
            "Ca",
            dimensionedScalar("one", inv(dimTime), 1.0)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::aggregationKernel::~aggregationKernel()
{}


Foam::scalar
Foam::populationBalanceSubModels::aggregationKernel::aggregationSource
(
    const label momentOrder,
    const label celli,
    const scalarQuadratureApproximation& quadrature,
    const label enviroment
)
{
    scalar aSource = 0.0;

    const PtrList<volNode>& nodes = quadrature.nodes();

    if (!nodes[0].extended())   // Non-extended quadrature case
    {
        forAll(nodes, pNode1i)
        {
            const volNode& node1 = nodes[pNode1i];
            const volScalarField& pWeight1 = node1.primaryWeight();
            const volScalarField& pAbscissa1 = node1.primaryAbscissae()[0];

            forAll(nodes, pNode2i)
            {
                const volNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();
                const volScalarField& pAbscissa2 = node2.primaryAbscissae()[0];

                // Remove small negative values in abscissae
                scalar bAbscissa1 = max(pAbscissa1[celli], 0.0);
                scalar bAbscissa2 = max(pAbscissa2[celli], 0.0);

                aSource +=
                    pWeight1[celli]*
                    (
                        pWeight2[celli]*
                        (
                            0.5*pow // Birth
                            (
                                pow3(bAbscissa1) + pow3(bAbscissa2),
                                momentOrder/3.0
                            )
                            - pow(bAbscissa1, momentOrder)
                        )
                       *Ka
                        (
                            bAbscissa1, bAbscissa2, celli, enviroment
                        )
                    );
            }
        }

        return aSource;
    }

    forAll(nodes, pNode1i)      // Extended quadrature case
    {
        const volNode& node1 = nodes[pNode1i];
        const volScalarField& pWeight1 = node1.primaryWeight();

        forAll(node1.secondaryWeights(), sNode1i)
        {
            const volScalarField& sWeight1 = node1.secondaryWeights()[sNode1i];

            const volScalarField& sAbscissa1
                = node1.secondaryAbscissae()[0][sNode1i];

            forAll(nodes, pNode2i)
            {
                const volNode& node2 = nodes[pNode2i];
                const volScalarField& pWeight2 = node2.primaryWeight();

                forAll(node2.secondaryWeights(), sNode2i)
                {
                    const volScalarField& sWeight2
                        = node2.secondaryWeights()[sNode2i];

                    const volScalarField& sAbscissa2
                        = node2.secondaryAbscissae()[0][sNode2i];

                    // Remove small negative values in abscissae
                    scalar bAbscissa1 = max(sAbscissa1[celli], 0.0);
                    scalar bAbscissa2 = max(sAbscissa2[celli], 0.0);

                    aSource +=
                        pWeight1[celli]*sWeight1[celli]*
                        (
                            pWeight2[celli]*sWeight2[celli]*
                            (
                                0.5*pow // Birth
                                (
                                    pow3(bAbscissa1) + pow3(bAbscissa2),
                                    momentOrder/3.0
                                )
                              - pow(bAbscissa1, momentOrder)
                            )
                           *Ka
                            (
                                bAbscissa1, bAbscissa2, celli, enviroment
                            )
                        );
                }
            }
        }
    }

    return aSource;
}

// ************************************************************************* //
