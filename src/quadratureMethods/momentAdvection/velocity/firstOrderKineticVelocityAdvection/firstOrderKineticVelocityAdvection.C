/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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

#include "firstOrderKineticVelocityAdvection.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace velocityAdvection
{
    defineTypeNameAndDebug(firstOrderKinetic, 0);

    addToRunTimeSelectionTable
    (
        velocityMomentAdvection,
        firstOrderKinetic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityAdvection::firstOrderKinetic::firstOrderKinetic
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const List<supportType>& supports
)
:
    velocityMomentAdvection(dict, quadrature, supports),
    weightScheme_("upwind"),
    scalarAbscissaeScheme_("upwind"),
    velocityAbscissaeScheme_("upwind")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityAdvection::firstOrderKinetic::~firstOrderKinetic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityAdvection::firstOrderKinetic::interpolateNodes()
{
    IStringStream weightOwnLimiter(weightScheme_);
    IStringStream scalarAbscissaeOwnLimiter(scalarAbscissaeScheme_);
    IStringStream velocityAbscissaeOwnLimiter(velocityAbscissaeScheme_);

    tmp<surfaceInterpolationScheme<scalar>> weightOwnScheme
    (
        fvc::scheme<scalar>(own_, weightOwnLimiter)
    );

    tmp<surfaceInterpolationScheme<scalar>> scalarAbscissaeOwnScheme
    (
        fvc::scheme<scalar>
        (
            own_,
            scalarAbscissaeOwnLimiter
        )
    );

    tmp<surfaceInterpolationScheme<vector>> velocityAbscissaeOwnScheme
    (
        fvc::scheme<vector>
        (
            own_,
            velocityAbscissaeOwnLimiter
        )
    );

    IStringStream weightNeiLimiter(weightScheme_);
    IStringStream scalarAbscissaeNeiLimiter(scalarAbscissaeScheme_);
    IStringStream velocityAbscissaeNeiLimiter(velocityAbscissaeScheme_);

    tmp<surfaceInterpolationScheme<scalar>> weightNeiScheme
    (
        fvc::scheme<scalar>(nei_, weightNeiLimiter)
    );

    tmp<surfaceInterpolationScheme<scalar>> scalarAbscissaeNeiScheme
    (
        fvc::scheme<scalar>
        (
            nei_,
            scalarAbscissaeNeiLimiter
        )
    );

    tmp<surfaceInterpolationScheme<vector>> velocityAbscissaeNeiScheme
    (
        fvc::scheme<vector>
        (
            nei_,
            velocityAbscissaeNeiLimiter
        )
    );

    PtrList<surfaceVelocityNode>& nodesNei = nodesNei_();
    PtrList<surfaceVelocityNode>& nodesOwn = nodesOwn_();

    forAll(nodes_, nodei)
    {
        const volVelocityNode& node(nodes_[nodei]);
        surfaceVelocityNode& nodeNei(nodesNei[nodei]);
        surfaceVelocityNode& nodeOwn(nodesOwn[nodei]);

        nodeOwn.weight() =
            weightOwnScheme().interpolate(node.weight());

        nodeOwn.velocityAbscissae() =
            velocityAbscissaeOwnScheme().interpolate(node.velocityAbscissae());

        nodeNei.weight() =
            weightNeiScheme().interpolate(node.weight());

        nodeNei.velocityAbscissae() =
            velocityAbscissaeNeiScheme().interpolate(node.velocityAbscissae());

        forAll(node.abscissae(), cmpt)
        {
            nodeOwn.abscissae()[cmpt] =
                scalarAbscissaeOwnScheme().interpolate
                (
                    node.abscissae()[cmpt]
                );

            nodeNei.abscissae()[cmpt] =
                scalarAbscissaeNeiScheme().interpolate
                (
                    node.abscissae()[cmpt]
                );
        }
    }
}


Foam::scalar
Foam::velocityAdvection::firstOrderKinetic::realizableCo() const
{
    const fvMesh& mesh = this->own_.mesh();

    scalarField maxCoNum(mesh.nCells(), scalar(1));

    forAll(this->nodes_, nodei)
    {
        surfaceScalarField phiOwn
        (
            mag(this->nodesOwn_()[nodei].velocityAbscissae() & mesh.Sf())
        );

        surfaceScalarField phiNei
        (
            mag(this->nodesNei_()[nodei].velocityAbscissae() & mesh.Sf())
        );

        forAll(moments_[0], celli)
        {
            const labelList& cell = mesh.cells()[celli];

            scalar den = 0;

            forAll(cell, facei)
            {
                if (cell[facei] < mesh.nInternalFaces())
                {
                    den +=
                        max
                        (
                            phiOwn[cell[facei]],
                            phiNei[cell[facei]]
                        );
                }
            }

            // The limit is taken once the sum of the fluxes leaving the
            // cell is complete. Taken inside the loop above it was still
            // right, the sum only growing, but the clamp beside it is an
            // assignment: a cell whose first faces contributed nothing
            // carried the floor in its denominator from there on.
            den = max(den, SMALL);

            maxCoNum[celli] =
                min
                (
                    maxCoNum[celli],
                    mesh.V()[celli]/(den*mesh.time().deltaTValue())
                );
        }
    }

    return gMin(maxCoNum);
}

Foam::scalar Foam::velocityAdvection::firstOrderKinetic::CoNum() const
{
    scalar CoNum = 0.0;
    const fvMesh& mesh = own_.mesh();
    forAll(nodes_, nodei)
    {
        CoNum =
            max
            (
                CoNum,
                0.5*gMax
                (
                    fvc::surfaceSum
                    (
                        mag(fvc::flux(nodes_[nodei].velocityAbscissae()))
                    )().primitiveField()/mesh.V().field()
                )*mesh.time().deltaTValue()
            );
    }

    return CoNum;
}

void Foam::velocityAdvection::firstOrderKinetic::resetMomentFluxes()
{
    forAll(divMoments_, divi)
    {
        divMoments_[divi] =
            dimensionedScalar
            (
                "0",
                moments_[divi].dimensions()/dimTime,
                Zero
            );
    }
}


void Foam::velocityAdvection::firstOrderKinetic::addNodeToMomentFluxes
(
    const label nodei,
    const surfaceScalarField& phiOwn,
    const surfaceScalarField& phiNei
)
{
    const dimensionedScalar zeroPhi("zero", dimVolume/dimTime, Zero);

    const labelList& scalarIndexes = nodes_[0].scalarIndexes();
    const labelList& velocityIndexes = nodes_[0].velocityIndexes();

    const surfaceVelocityNode& nodeNei(nodesNei_()[nodei]);
    const surfaceVelocityNode& nodeOwn(nodesOwn_()[nodei]);

    const surfaceScalarField& weightOwn = nodeOwn.weight();
    const surfaceScalarField& weightNei = nodeNei.weight();

    const PtrList<surfaceScalarField>& scalarAbscissaeOwn =
        nodeOwn.abscissae();

    const PtrList<surfaceScalarField>& scalarAbscissaeNei =
        nodeNei.abscissae();

    const surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
    const surfaceVectorField& UNei = nodeNei.velocityAbscissae();

    forAll(divMoments_, divi)
    {
        const labelList& momentOrder = momentOrders_[divi];

        surfaceScalarField momentCmptOwn(weightOwn);
        surfaceScalarField momentCmptNei(weightNei);

        // The moment is the weight times the abscissa of each coordinate
        // raised to its own order. The dimensions of the accumulator are
        // reset on the way because that order varies from one moment to
        // the next, which the dimension check of a field cannot express.
        forAll(scalarIndexes, cmpti)
        {
            const label cmpt = scalarIndexes[cmpti];
            const label cmptMomentOrder = momentOrder[cmpt];

            if (cmptMomentOrder > 0)
            {
                const surfaceScalarField& abscissaOwnCmpt =
                    scalarAbscissaeOwn[cmpti];

                const surfaceScalarField& abscissaNeiCmpt =
                    scalarAbscissaeNei[cmpti];

                tmp<surfaceScalarField> mOwnPow =
                    momentCmptOwn*pow(abscissaOwnCmpt, cmptMomentOrder);

                tmp<surfaceScalarField> mNeiPow =
                    momentCmptNei*pow(abscissaNeiCmpt, cmptMomentOrder);

                momentCmptOwn.dimensions().reset(mOwnPow().dimensions());
                momentCmptOwn == mOwnPow;

                momentCmptNei.dimensions().reset(mNeiPow().dimensions());
                momentCmptNei == mNeiPow;
            }
        }

        forAll(velocityIndexes, cmpti)
        {
            const label cmpt = velocityIndexes[cmpti];
            const label cmptMomentOrder = momentOrder[cmpt];

            if (cmptMomentOrder > 0)
            {
                tmp<surfaceScalarField> abscissaOwnCmpt =
                    UOwn.component(cmpti);

                tmp<surfaceScalarField> abscissaNeiCmpt =
                    UNei.component(cmpti);

                tmp<surfaceScalarField> mOwnPow =
                    momentCmptOwn*pow(abscissaOwnCmpt, cmptMomentOrder);

                tmp<surfaceScalarField> mNeiPow =
                    momentCmptNei*pow(abscissaNeiCmpt, cmptMomentOrder);

                momentCmptOwn.dimensions().reset(mOwnPow().dimensions());
                momentCmptOwn == mOwnPow;

                momentCmptNei.dimensions().reset(mNeiPow().dimensions());
                momentCmptNei == mNeiPow;
            }
        }

        divMoments_[divi] +=
            fvc::surfaceIntegrate
            (
                momentCmptOwn*max(phiOwn, zeroPhi)
              + momentCmptNei*min(phiNei, zeroPhi)
            );
    }
}


void Foam::velocityAdvection::firstOrderKinetic::update()
{
    const fvMesh& mesh = own_.mesh();

    // Interpolate weights and abscissae
    interpolateNodes();

    // Set velocities at boundaries for rebounding
    updateBoundaryConditions();

    resetMomentFluxes();

    forAll(nodes_, nodei)
    {
        // The flux is carried by the velocity abscissa of the node itself
        surfaceScalarField phiOwn
        (
            nodesOwn_()[nodei].velocityAbscissae() & mesh.Sf()
        );

        surfaceScalarField phiNei
        (
            nodesNei_()[nodei].velocityAbscissae() & mesh.Sf()
        );

        addNodeToMomentFluxes(nodei, phiOwn, phiNei);
    }
}


void Foam::velocityAdvection::firstOrderKinetic::update
(
    const surfaceScalarField& phi,
    const bool wallCollisions
)
{
    // Interpolate weights and abscissae
    interpolateNodes();

    // Set velocities at boundaries for rebounding
    if (wallCollisions)
    {
        updateBoundaryConditions();
    }

    resetMomentFluxes();

    forAll(nodes_, nodei)
    {
        // Every node is carried by the flux that was given, so the two
        // sides of a face see the same one
        addNodeToMomentFluxes(nodei, phi, phi);
    }
}


void Foam::velocityAdvection::firstOrderKinetic::update
(
    const mappedPtrList<volVectorField>& Us,
    const bool wallCollisions
)
{
    const fvMesh& mesh = own_.mesh();

    // Interpolate weights and abscissae
    interpolateNodes();

    IStringStream velocityAbscissaeOwnLimiter(velocityAbscissaeScheme_);

    tmp<surfaceInterpolationScheme<vector>> velocityAbscissaeOwnScheme
    (
        fvc::scheme<vector>
        (
            own_,
            velocityAbscissaeOwnLimiter
        )
    );

    IStringStream velocityAbscissaeNeiLimiter(velocityAbscissaeScheme_);

    tmp<surfaceInterpolationScheme<vector>> velocityAbscissaeNeiScheme
    (
        fvc::scheme<vector>
        (
            nei_,
            velocityAbscissaeNeiLimiter
        )
    );

    // Set velocities at boundaries for rebounding
    if (wallCollisions)
    {
        updateBoundaryConditions();
    }

    resetMomentFluxes();

    forAll(nodes_, nodei)
    {
        // The flux is carried by the velocity field given for the node
        // rather than by its own abscissa
        surfaceVectorField VOwn
        (
            velocityAbscissaeOwnScheme().interpolate(Us[nodei])
        );

        surfaceVectorField VNei
        (
            velocityAbscissaeNeiScheme().interpolate(Us[nodei])
        );

        surfaceScalarField phiOwn(VOwn & mesh.Sf());
        surfaceScalarField phiNei(VNei & mesh.Sf());

        addNodeToMomentFluxes(nodei, phiOwn, phiNei);
    }
}


// ************************************************************************* //
