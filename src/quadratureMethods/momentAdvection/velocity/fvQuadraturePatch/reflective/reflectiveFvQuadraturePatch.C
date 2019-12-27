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

#include "reflectiveFvQuadraturePatch.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reflectiveFvQuadraturePatch, 0);

    addToRunTimeSelectionTable
    (
        fvQuadraturePatch,
        reflectiveFvQuadraturePatch,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reflectiveFvQuadraturePatch::reflectiveFvQuadraturePatch
(
    const fvPatch& patch,
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    PtrList<surfaceVelocityNode>& nodesOwn,
    PtrList<surfaceVelocityNode>& nodesNei
)
:
    fvQuadraturePatch(patch, dict, quadrature, nodesOwn, nodesNei),
    ew_(readScalar(dict.lookup("e")))
{
    if (!isA<wallFvPatch>(patch_))
    {
        FatalErrorInFunction
            << "Wall physical boundary required, but type "
            << patch_.type() << " specified."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reflectiveFvQuadraturePatch::~reflectiveFvQuadraturePatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::reflectiveFvQuadraturePatch::wallTangentVelocity
(
    const vectorField& U,
    const vectorField& n
) const
{
    return (U - n*(U & n));
}


void Foam::reflectiveFvQuadraturePatch::update()
{
    if (!patch_.size())
    {
        return;
    }

    const vectorField& bfSf(patch_.Sf());
    vectorField bfNorm(patch_.nf());

    scalarField Gin(bfSf.size(), Zero);
    scalarField Gout(bfSf.size(), Zero);

    forAll(quadrature_.nodes(), nodei)
    {
        const volVelocityNode& node = quadrature_.nodes()[nodei];
        surfaceVelocityNode& nodeNei(nodesNei_[nodei]);
        surfaceVelocityNode& nodeOwn(nodesOwn_[nodei]);

        const volScalarField& weight = node.primaryWeight();
        surfaceScalarField& weightOwn = nodeOwn.primaryWeight();
        surfaceScalarField& weightNei = nodeNei.primaryWeight();
        const volVectorField& U = node.velocityAbscissae();
        surfaceVectorField& UOwn = nodeOwn.velocityAbscissae();
        surfaceVectorField& UNei = nodeNei.velocityAbscissae();

        scalarField& bfwOwn = weightOwn.boundaryFieldRef()[patchi_];
        scalarField& bfwNei = weightNei.boundaryFieldRef()[patchi_];
        vectorField& bfUOwn = UOwn.boundaryFieldRef()[patchi_];
        vectorField& bfUNei = UNei.boundaryFieldRef()[patchi_];

        bfwOwn = weight.boundaryField()[patchi_].patchInternalField();
        bfwNei = bfwOwn;

        bfUOwn = U.boundaryField()[patchi_].patchInternalField();

        tmp<vectorField> vn
        (
            -this->ew_*(bfUOwn & bfNorm)*bfNorm
        );

        bfUNei = vn + wallTangentVelocity(bfUOwn, bfNorm);

        Gin += max(scalar(0), bfUOwn & bfSf)*bfwOwn;
        Gout -= min(scalar(0), bfUNei & bfSf)*bfwNei;
    }

    //- Scale to ensure zero flux
    if (this->ew_ < 1)
    {
        scalarField weightScale(Gin/(Gout + SMALL));

        forAll(quadrature_.nodes(), nodei)
        {
            scalarField& bfWNei =
                nodesNei_[nodei].primaryWeight().boundaryFieldRef()[patchi_];

            bfWNei *= weightScale;
        }
    }
}


// ************************************************************************* //
