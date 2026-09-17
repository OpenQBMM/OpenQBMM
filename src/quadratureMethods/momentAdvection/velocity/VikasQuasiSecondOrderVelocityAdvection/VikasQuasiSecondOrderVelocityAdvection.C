/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2025 Alberto Passalacqua
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
    const List<supportType>& supports
)
:
    firstOrderKinetic(dict, quadrature, supports),
    minWeight_(dict.lookupOrDefault<scalar>("minWeight", 1.0e-6)),
    nDroppedNodes_(0),
    lowestDroppedLimit_(GREAT),
    warnedDroppedNodes_(false)
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

            scalar num = this->nodes_[nodei].weight()[celli];
            scalar den = 0;
            forAll(cell, facei)
            {
                if (cell[facei] < mesh.nInternalFaces())
                {
                    if (own[cell[facei]] == celli)
                    {
                        den +=
                            this->nodesOwn_()[nodei].weight()[cell[facei]]
                           *max(phiOwn[cell[facei]], scalar(0));
                    }
                    else if (nei[cell[facei]] == celli)
                    {
                        den -=
                            this->nodesNei_()[nodei].weight()[cell[facei]]
                           *min(phiNei[cell[facei]], scalar(0));
                    }
                }
            }

            // The flux leaving through the boundary, with the weight of
            // the node on the inner side of those faces
            {
                const surfaceScalarField::Boundary& phiBf =
                    phiOwn.boundaryField();

                const surfaceScalarField::Boundary& wBf =
                    this->nodesOwn_()[nodei].weight().boundaryField();

                forAll(cell, facei)
                {
                    if (cell[facei] >= mesh.nInternalFaces())
                    {
                        const label patchi =
                            mesh.boundaryMesh().whichPatch(cell[facei]);

                        if (patchi < 0)
                        {
                            continue;
                        }

                        const label pFacei =
                            cell[facei] - mesh.boundaryMesh()[patchi].start();

                        if (pFacei < phiBf[patchi].size())
                        {
                            den +=
                                wBf[patchi][pFacei]
                               *max(phiBf[patchi][pFacei], scalar(0));
                        }
                    }
                }
            }

            // As in the scheme this derives from, the limit is taken
            // once the sum of the fluxes leaving the cell is complete.
            // Taken inside the loop above it was still right, the sum only
            // growing, but the clamp beside it is an assignment: a cell
            // whose first faces contributed nothing carried the floor in
            // its denominator from there on.
            den = max(den, SMALL);

            const scalar limit =
                num*mesh.V()[celli]/(den*mesh.time().deltaTValue());

            if (num > minWeight_)
            {
                maxCoNum[celli] = min(maxCoNum[celli], limit);
            }
            else if
            (
                limit < maxCoNum[celli]
             && den*mesh.time().deltaTValue() > minWeight_*mesh.V()[celli]
            )
            {
                // A node too light to be counted would have lowered the
                // limit, and what it carries out over the step is more
                // than the weight it was dropped for: it loses more than
                // it has. There is nothing to do for it in this scheme, so
                // it is counted and reported once. A node whose outflow is
                // itself below that weight is empty and left alone.
                nDroppedNodes_++;
                lowestDroppedLimit_ = min(lowestDroppedLimit_, limit);
            }
        }
    }

    if (nDroppedNodes_ > 0 && !warnedDroppedNodes_)
    {
        WarningInFunction
            << "Nodes with a weight below minWeight " << minWeight_
            << " were left out of the Courant limit while they would have"
            << " lowered it, in " << nDroppedNodes_
            << " instances so far, to as little as " << lowestDroppedLimit_
            << " of the step." << nl
            << "    This scheme has no treatment for them; the limit it"
            << " reports does not cover what they carry." << endl;

        warnedDroppedNodes_ = true;
    }

    return gMin(maxCoNum);
}

// ************************************************************************* //
