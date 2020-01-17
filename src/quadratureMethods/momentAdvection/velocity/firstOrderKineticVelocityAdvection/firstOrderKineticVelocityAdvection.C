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
    Copyright (C) 2019 Alberto Passalacqua
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
    const word& support
)
:
    velocityMomentAdvection(dict, quadrature, support),
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

        nodeOwn.primaryWeight() =
            weightOwnScheme().interpolate(node.primaryWeight());

        nodeOwn.velocityAbscissae() =
            velocityAbscissaeOwnScheme().interpolate(node.velocityAbscissae());

        nodeNei.primaryWeight() =
            weightNeiScheme().interpolate(node.primaryWeight());

        nodeNei.velocityAbscissae() =
            velocityAbscissaeNeiScheme().interpolate(node.velocityAbscissae());

        forAll(node.primaryAbscissae(), cmpt)
        {
            nodeOwn.primaryAbscissae()[cmpt] =
                scalarAbscissaeOwnScheme().interpolate
                (
                    node.primaryAbscissae()[cmpt]
                );

            nodeNei.primaryAbscissae()[cmpt] =
                scalarAbscissaeNeiScheme().interpolate
                (
                    node.primaryAbscissae()[cmpt]
                );
        }
    }
}


Foam::scalar
Foam::velocityAdvection::firstOrderKinetic::realizableCo() const
{
    const fvMesh& mesh = this->own_.mesh();
    surfaceVectorField Sf(mesh.Sf());

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

                den = max(den, SMALL);

                maxCoNum[celli] =
                    min
                    (
                        maxCoNum[celli],
                        mesh.V()[celli]
                       /(den*mesh.time().deltaTValue())
                    );
            }
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

void Foam::velocityAdvection::firstOrderKinetic::update()
{
    const fvMesh& mesh = own_.mesh();
    dimensionedScalar zeroPhi("zero", dimVolume/dimTime, Zero);

    // Interpolate weights and abscissae
    interpolateNodes();

    // Set velocities at boundaries for rebounding
    updateBoundaryConditions();

    // Zero moment flux
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

    const labelList& scalarIndexes = nodes_[0].scalarIndexes();
    const labelList& velocityIndexes = nodes_[0].velocityIndexes();

    forAll(nodes_, nodei)
    {
        const surfaceVelocityNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVelocityNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();

        const PtrList<surfaceScalarField>& scalarAbscissaeOwn =
            nodeOwn.primaryAbscissae();

        const PtrList<surfaceScalarField>& scalarAbscissaeNei =
            nodeNei.primaryAbscissae();

        const surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
        const surfaceVectorField& UNei = nodeNei.velocityAbscissae();

        surfaceScalarField phiOwn(UOwn & mesh.Sf());
        surfaceScalarField phiNei(UNei & mesh.Sf());

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];

            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

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
                        momentCmptOwn
                       *pow
                        (
                            abscissaOwnCmpt,
                            cmptMomentOrder
                        );

                    tmp<surfaceScalarField> mNeiPow =
                        momentCmptNei
                       *pow
                        (
                            abscissaNeiCmpt,
                            cmptMomentOrder
                        );

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
                        momentCmptOwn
                       *pow
                        (
                            abscissaOwnCmpt,
                            cmptMomentOrder
                        );

                    tmp<surfaceScalarField> mNeiPow =
                        momentCmptNei
                       *pow
                        (
                            abscissaNeiCmpt,
                            cmptMomentOrder
                        );

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
}

void Foam::velocityAdvection::firstOrderKinetic::update
(
    const surfaceScalarField& phi,
    const bool wallCollisions
)
{
    dimensionedScalar zeroPhi("zero", dimVolume/dimTime, Zero);

    // Interpolate weights and abscissae
    interpolateNodes();

    // Set velocities at boundaries for rebounding
    if (wallCollisions)
    {
        updateBoundaryConditions();
    }

    // Zero moment fluxes
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

    const labelList& scalarIndexes = nodes_[0].scalarIndexes();
    const labelList& velocityIndexes = nodes_[0].velocityIndexes();

    forAll(nodes_, nodei)
    {
        const surfaceVelocityNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVelocityNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();

        const PtrList<surfaceScalarField>& scalarAbscissaeOwn =
            nodeOwn.primaryAbscissae();

        const PtrList<surfaceScalarField>& scalarAbscissaeNei =
            nodeNei.primaryAbscissae();

        const surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
        const surfaceVectorField& UNei = nodeNei.velocityAbscissae();

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];

            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

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
                        momentCmptOwn
                       *pow
                        (
                            abscissaOwnCmpt,
                            cmptMomentOrder
                        );

                    tmp<surfaceScalarField> mNeiPow =
                        momentCmptNei
                       *pow
                        (
                            abscissaNeiCmpt,
                            cmptMomentOrder
                        );

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
                        momentCmptOwn
                       *pow
                        (
                            abscissaOwnCmpt,
                            cmptMomentOrder
                        );

                    tmp<surfaceScalarField> mNeiPow =
                        momentCmptNei
                       *pow
                        (
                            abscissaNeiCmpt,
                            cmptMomentOrder
                        );

                    momentCmptOwn.dimensions().reset(mOwnPow().dimensions());
                    momentCmptOwn == mOwnPow;

                    momentCmptNei.dimensions().reset(mNeiPow().dimensions());
                    momentCmptNei == mNeiPow;
                }
            }

            divMoments_[divi] +=
                fvc::surfaceIntegrate
                (
                    momentCmptOwn*max(phi, zeroPhi)
                  + momentCmptNei*min(phi, zeroPhi)
                );
        }
    }
}

void Foam::velocityAdvection::firstOrderKinetic::update
(
    const mappedPtrList<volVectorField>& Us,
    const bool wallCollisions
)
{
    const fvMesh& mesh = own_.mesh();
    dimensionedScalar zeroPhi("zero", dimVolume/dimTime, Zero);

    // Interplate weights and abscissae
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

    // Zero moment fluxes
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

    const labelList& scalarIndexes = nodes_[0].scalarIndexes();
    const labelList& velocityIndexes = nodes_[0].velocityIndexes();

    forAll(nodes_, nodei)
    {
        const surfaceVelocityNode& nodeNei(nodesNei_()[nodei]);
        const surfaceVelocityNode& nodeOwn(nodesOwn_()[nodei]);

        const surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        const surfaceScalarField& weightNei = nodeNei.primaryWeight();

        const PtrList<surfaceScalarField>& scalarAbscissaeOwn =
            nodeOwn.primaryAbscissae();

        const PtrList<surfaceScalarField>& scalarAbscissaeNei =
            nodeNei.primaryAbscissae();

        const surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
        const surfaceVectorField& UNei = nodeNei.velocityAbscissae();

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

        forAll(divMoments_, divi)
        {
            const labelList& momentOrder = momentOrders_[divi];

            surfaceScalarField momentCmptOwn(weightOwn);
            surfaceScalarField momentCmptNei(weightNei);

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
                        momentCmptOwn
                       *pow
                        (
                            abscissaOwnCmpt,
                            cmptMomentOrder
                        );

                    tmp<surfaceScalarField> mNeiPow =
                        momentCmptNei
                       *pow
                        (
                            abscissaNeiCmpt,
                            cmptMomentOrder
                        );

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
                        momentCmptOwn
                       *pow
                        (
                            abscissaOwnCmpt,
                            cmptMomentOrder
                        );

                    tmp<surfaceScalarField> mNeiPow =
                        momentCmptNei
                       *pow
                        (
                            abscissaNeiCmpt,
                            cmptMomentOrder
                        );
                        
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
}

// ************************************************************************* //
