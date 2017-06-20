/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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

#include "zeta.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zeta, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentAdvection,
        zeta,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeta::zeta
(
    const dictionary& dict,
    const univariateQuadratureApproximation& quadrature,
    const surfaceScalarField& phi,
    const word& support
)
:
    univariateMomentAdvection(dict, quadrature, phi, support),
    nZetas_(nMoments_ - 1),
    zetas_(nZetas_),
    zetasNei_(nZetas_),
    zetasOwn_(nZetas_),
    momentsNei_(nMoments_),
    momentsOwn_(nMoments_),
    phi_(phi)
{
    // Populating zeta_k fields and interpolated zeta_k fields
    forAll(zetas_, zetai)
    {
        zetas_.set
        (
            zetai,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName(name_, "zeta") + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        zetasNei_.set
        (
            zetai,
            new surfaceScalarField
            (
                IOobject
                (
                    "zetaNei" + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        zetasOwn_.set
        (
            zetai,
            new surfaceScalarField
            (
                IOobject
                (
                    "zetaOwn" + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }

    // Setting face values of moments
    forAll(momentsNei_, momenti)
    {
        momentsNei_.set
        (
            momenti,
            new surfaceScalarField
            (
                "momentNei" + Foam::name(momenti),
                fvc::interpolate(moments_[momenti])
            )
        );

        momentsOwn_.set
        (
            momenti,
            new surfaceScalarField
            (
                "momentOwn" + Foam::name(momenti),
                fvc::interpolate(moments_[momenti])
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zeta::~zeta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeta::interpolateZetas()
{
    forAll(zetas_, zetai)
    {
        zetasNei_[zetai] =
            fvc::interpolate(zetas_[zetai], nei_, "reconstruct(zeta)");

        zetasOwn_[zetai] =
            fvc::interpolate(zetas_[zetai], own_, "reconstruct(zeta)");
    }
}

void Foam::zeta::zetaToMoments
(
    const scalarList& zetaf,
    scalarList& mf
)
{
    scalarSquareMatrix S(nMoments_, 0.0);

    for (label i = 0; i < nZetas_; i++)
    {
        S[0][i] = 1.0;
    }

    for (label i = 1; i < nZetas_; i++)
    {
        for (label j = i; j < nZetas_; j++)
        {
            S[i][j] = S[i][j - 1] + zetaf[j - i]*S[i - 1][j];
        }
    }

    scalarList prod(nMoments_, 1.0);

    prod[1] = zetaf[0];

    for (label i = 2; i < nZetas_; i++)
    {
        prod[i] = prod[i - 1]*zetaf[i - 1];
    }

    // Resetting moments to zero
    mf = 0.0;

    // Computing moments
    mf[0] = 1.0;
    mf[1] = zetaf[0];

    for (label i = 2; i < nMoments_; i++)
    {
        for (label j = 0; j <= i/2; j++)
        {
            mf[i] += prod[i - 2*j]*sqr(S[j][i - j]);
        }
    }
}

void Foam::zeta::computeZetaFields()
{
    // Cell-center values
    forAll(moments_[0], celli)
    {
        if (moments_[0][celli] >= SMALL)
        {
            univariateMomentSet m(nMoments_, support_);

            for (label mi = 0; mi < nMoments_; mi++)
            {
                m[mi] = moments_[mi][celli];
            }

            scalarList zetas(m.zetas());

            for (label zetai = 0; zetai < nZetas_; zetai++)
            {
                zetas_[zetai][celli] = zetas[zetai];
            }
        }
    }

    // Boundary conditions
    const volScalarField::Boundary& bf = zetas_[0].boundaryField();

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            if (moments_[0].boundaryField()[patchi][facei] >= SMALL)
            {
                univariateMomentSet m(nMoments_, support_);

                for (label mi = 0; mi < nMoments_; mi++)
                {
                    m[mi] = moments_[mi].boundaryField()[patchi][facei];
                }

                scalarList zetas(m.zetas());

                for (label zetai = 0; zetai < nZetas_; zetai++)
                {
                    volScalarField& zi = zetas_[zetai];
                    volScalarField::Boundary& ziBf = zi.boundaryFieldRef();
                    ziBf[patchi][facei] = zetas[zetai];
                }
            }
        }
    }
}

Foam::scalar Foam::zeta::realizableCo()
{
    const fvMesh& mesh(phi_.mesh());
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    scalarField internalCo(mesh.nCells(), 0.0);

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        if (phi_[facei] > 0)
        {
            internalCo[own[facei]] += 1;
        }
        else if (phi_[facei] < 0)
        {
            internalCo[nei[facei]] += 1;
        }
    }

    internalCo = 1.0/(internalCo + 1.0);

    Info << internalCo << endl;

    return min(gMin(internalCo), 1.0/3.0);
}

void Foam::zeta::update()
{
    // Reconstructing zero-order moment on cell faces (standard min-mod).
    surfaceScalarField m0Own
    (
        "m0Own",
        fvc::interpolate(moments_[0], own_, "reconstruct(m0)")
    );

    surfaceScalarField m0Nei
    (
        "m0Nei",
        fvc::interpolate(moments_[0], nei_, "reconstruct(m0)")
    );

    // Compute zeta fields
    computeZetaFields();

    // Reconstructing zeta_k on cell faces
    interpolateZetas();

    // Recompute moments at sides of cell faces
    updateMomentFieldsFromZetas(m0Nei, zetasNei_, momentsNei_);
    updateMomentFieldsFromZetas(m0Own, zetasOwn_, momentsOwn_);

    // Calculate moment advection term
    dimensionedScalar zeroPhi("zero", phi_.dimensions(), 0.0);

    forAll(divMoments_, divi)
    {
        volScalarField divMoment
        (
            IOobject
            (
                "divMoment",
                moments_[0].mesh().time().timeName(),
                moments_[0].mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            moments_[0].mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        );

        surfaceScalarField mFlux
        (
            momentsNei_[divi]*min(phi_, zeroPhi)
          + momentsOwn_[divi]*max(phi_, zeroPhi)
        );

        fvc::surfaceIntegrate(divMoment.ref(), mFlux);
        divMoment.ref().dimensions().reset(moments_[divi].dimensions()/dimTime);

        divMoments_[divi].replace(0, divMoment);
    }
}

void Foam::zeta::updateMomentFieldsFromZetas
(
    const surfaceScalarField m0f,
    const PtrList<surfaceScalarField>& zetaf,
    PtrList<surfaceScalarField>& mf
)
{
    forAll(zetaf[0], facei)
    {
        scalarList zf(nZetas_);

        for (label zetai = 0; zetai < nZetas_; zetai++)
        {
            zf[zetai] = zetaf[zetai][facei];
        }

        scalarList mFace(nMoments_, 0.0);
        zetaToMoments(zf, mFace);

        for (label mi = 0; mi < nMoments_; mi++)
        {
            mf[mi][facei] = m0f[facei]*mFace[mi];
        }
    }

    // Boundary conditions
    const surfaceScalarField::Boundary& bf = zetaf[0].boundaryField();

    forAll(bf, patchi)
    {
        const fvsPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            scalarList zf(nZetas_);

            for (label zetai = 0; zetai < nZetas_; zetai++)
            {
                zf[zetai] = zetaf[zetai].boundaryField()[patchi][facei];
            }

            scalarList mFace(nMoments_, 0.0);
            zetaToMoments(zf, mFace);

            for (label mi = 0; mi < nMoments_; mi++)
            {
                mf[mi].boundaryFieldRef()[patchi][facei] = m0f.boundaryField()[patchi][facei]*mFace[mi];
            }
        }
    }
}

// ************************************************************************* //
