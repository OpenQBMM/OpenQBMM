/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 Alberto Passalacqua
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

#include "VikasQuasiSecondOrderVelocityAdvection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{
    defineTypeNameAndDebug(VikasQuasiSecondOrder, 0);

    addToRunTimeSelectionTable
    (
        velocityMomentAdvection,
        VikasQuasiSecondOrder,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityAdvection::VikasQuasiSecondOrder::VikasQuasiSecondOrder
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const word& support
)
:
    firstOrderKinetic(dict, quadrature, support)
{
    weightScheme_ = "Minmod";
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityAdvection::VikasQuasiSecondOrder::~VikasQuasiSecondOrder()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::velocityAdvection::VikasQuasiSecondOrder::realizableCo() const
{
    const fvMesh& mesh = this->own_.mesh();
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    surfaceVectorField Sf(mesh.Sf());

    scalarField maxCoNum(mesh.nCells(), scalar(1));

    forAll(this->nodes_, nodei)
    {

        surfaceScalarField phiOwn
        (
            this->nodesOwn_()[nodei].velocityAbscissae() & mesh.Sf()
        );

        surfaceScalarField phiNei
        (
            this->nodesNei_()[nodei].velocityAbscissae() & mesh.Sf()
        );

        forAll(moments_[0], celli)
        {
            const labelList& cell = mesh.cells()[celli];

            scalar num = this->nodes_[nodei].primaryWeight()[celli];
            scalar den = 0;
            forAll(cell, facei)
            {
                if (cell[facei] < mesh.nInternalFaces())
                {
                    if (own[cell[facei]] == celli)
                    {
                        den +=
                            this->nodesOwn_()[nodei].primaryWeight()[cell[facei]]
                           *max(phiOwn[cell[facei]], scalar(0));
                    }
                    else if (nei[cell[facei]] == celli)
                    {
                        den -=
                            this->nodesNei_()[nodei].primaryWeight()[cell[facei]]
                           *min(phiNei[cell[facei]], scalar(0));
                    }
                }
                if (num > 1e-6)
                {
                    den = max(den, SMALL);
                    maxCoNum[celli] =
                        min
                        (
                            maxCoNum[celli],
                            num*mesh.V()[celli]
                           /(den*mesh.time().deltaTValue())
                        );
                }
            }
        }
    }
    return gMin(maxCoNum);
}

// ************************************************************************* //
