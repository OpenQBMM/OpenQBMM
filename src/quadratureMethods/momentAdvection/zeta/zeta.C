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
#include "upwind.H"
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
    m0_(moments_[0]),
    m0Own_
    (
        "m0Own",
        fvc::interpolate(m0_, own_, "reconstruct(m0)")
    ),
    m0Nei_
    (
        "m0Nei",
        fvc::interpolate(m0_, nei_, "reconstruct(m0)")
    ),
    nZetas_(nMoments_ - 1),
    zetas_(nZetas_),
    zetasNei_(nZetas_),
    zetasOwn_(nZetas_),
    zetasUpwindNei_(nZetas_),
    zetasUpwindOwn_(nZetas_),
    zetasCorrNei_(nZetas_),
    zetasCorrOwn_(nZetas_),
    momentsNei_(nMoments_),
    momentsOwn_(nMoments_),
    nFacesOutgoingFlux_(m0_.mesh().nCells(), 0),
    nRealizableMoments_(m0_.mesh().nCells(), 0),
    nRealizableMomentsStar_(m0_.mesh().nCells(), 0),
    limiters_(nZetas_),
    cellLimiters_(nZetas_),
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
                    IOobject::groupName
                    (
                        IOobject::groupName("zeta", Foam::name(zetai)),
                        name_
                    ),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
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

        zetasUpwindNei_.set
        (
            zetai,
            new surfaceScalarField
            (
                IOobject
                (
                    "zetaUpwindNei" + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        zetasUpwindOwn_.set
        (
            zetai,
            new surfaceScalarField
            (
                IOobject
                (
                    "zetaUpwindOwn" + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        zetasCorrNei_.set
        (
            zetai,
            new surfaceScalarField
            (
                IOobject
                (
                    "zetaCorrNei" + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        zetasCorrOwn_.set
        (
            zetai,
            new surfaceScalarField
            (
                IOobject
                (
                    "zetaCorrOwn" + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        limiters_.set
        (
            zetai,
            new surfaceScalarField
            (
                IOobject
                (
                    "zetaLimiters" + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 1.0)
            )
        );

        cellLimiters_.set
        (
            zetai,
            new volScalarField
            (
                IOobject
                (
                    "zetaCellLimiters" + Foam::name(zetai),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, 1.0)
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

void Foam::zeta::interpolateFields()
{
    m0Own_ = fvc::interpolate(moments_[0], own_, "reconstruct(m0)");
    m0Nei_ = fvc::interpolate(moments_[0], nei_, "reconstruct(m0)");

    forAll(zetas_, zetai)
    {
        zetasNei_[zetai] =
            fvc::interpolate(zetas_[zetai], nei_, "reconstruct(zeta)");

        zetasOwn_[zetai] =
            fvc::interpolate(zetas_[zetai], own_, "reconstruct(zeta)");

        zetasUpwindNei_[zetai] =
            upwind<scalar>(zetas_[zetai].mesh(), nei_).flux(zetas_[zetai]);

        zetasUpwindOwn_[zetai] =
            upwind<scalar>(zetas_[zetai].mesh(), own_).flux(zetas_[zetai]);

        zetasCorrNei_[zetai] = zetasNei_[zetai] - zetasUpwindNei_[zetai];
        zetasCorrOwn_[zetai] = zetasOwn_[zetai] - zetasUpwindOwn_[zetai];
    }
}

void Foam::zeta::zetaToMoments
(
    const scalarList& zetaf,
    scalarList& mf,
    scalar m0
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

    if (m0 != 1.0)
    {
        for (label mi = 0; mi < nMoments_; mi++)
        {
            mf[mi] *= m0;
        }
    }
}

void Foam::zeta::computeZetaFields()
{
    // Cell-center values
    forAll(m0_, celli)
    {
        if (m0_[celli] >= SMALL)
        {
            univariateMomentSet m(nMoments_, support_);

            for (label mi = 0; mi < nMoments_; mi++)
            {
                m[mi] = moments_[mi][celli];
            }

            nRealizableMoments_[celli] = m.nRealizableMoments();

            scalarList zetas(m.zetas());

            for (label zetai = 0; zetai < nZetas_; zetai++)
            {
                zetas_[zetai][celli] = zetas[zetai];

                if (zetas_[zetai][celli] > 1.0e-7)
                {
                    zetas_[zetai][celli] = zetas[zetai];
                }
                else
                {
                    zetas_[zetai][celli] = 0.0;
                }
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
            if (m0_.boundaryField()[patchi][facei] >= SMALL)
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

void Foam::zeta::countFacesWithOutgoingFlux()
{
    const fvMesh& mesh(phi_.mesh());
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    nFacesOutgoingFlux_ = 0;

    // Counting internal faces with outgoing flux
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        if (phi_[facei] > 0)
        {
            nFacesOutgoingFlux_[own[facei]] += 1;
        }
        else if (phi_[facei] < 0)
        {
            nFacesOutgoingFlux_[nei[facei]] += 1;
        }
    }

    // Adding boundary faces with outgoing flux
    const surfaceScalarField::Boundary& phiBf = phi_.boundaryField();

    forAll(phiBf, patchi)
    {
        const fvsPatchScalarField& phiPf = phiBf[patchi];
        const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(phiPf, pFacei)
        {
            if (phiPf[pFacei] > 0)
            {
                nFacesOutgoingFlux_[pFaceCells[pFacei]] += 1;
            }
        }
    }
}

void Foam::zeta::limitZetas()
{
    const labelUList& owner = phi_.mesh().owner();
    const labelUList& neighb = phi_.mesh().neighbour();
    const scalarField& phiIf = phi_;
    const surfaceScalarField::Boundary& phiBf = phi_.boundaryField();
    const label nCells = phi_.mesh().nCells();
    const label nInternalFaces = phi_.mesh().nInternalFaces();

    countFacesWithOutgoingFlux();

    forAll(cellLimiters_, li)
    {
        forAll(cellLimiters_[0], celli)
        {
            cellLimiters_[li][celli] = 1.0;
        }
    }

    // First check on m* to identify cells in need of additional limitation
    scalarRectangularMatrix mPluses(nMoments_, nCells, 0.0);

    // Find m+ (moments reconstructed on cell faces with outgoing flux)
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = owner[facei];
        const label nei = neighb[facei];

        if (phi_[facei] > 0.0)
        {
            for (label mi = 0; mi < nMoments_; mi++)
            {
                mPluses[mi][own] += momentsOwn_[mi][facei];
            }
        }
        else
        {
            for (label mi = 0; mi < nMoments_; mi++)
            {
                mPluses[mi][nei] += momentsNei_[mi][facei];
            }
        }
    }

    // Adding boundary faces with outgoing flux
    forAll(phiBf, patchi)
    {
        const fvsPatchScalarField& phiPf = phiBf[patchi];

        const labelList& pFaceCells
            = phi_.mesh().boundary()[patchi].faceCells();

        forAll(phiPf, pFacei)
        {
            if (phiPf[pFacei] > 0)
            {
                for (label mi = 0; mi < nMoments_; mi++)
                {
                    mPluses[mi][pFaceCells[pFacei]] += momentsOwn_[mi][pFacei];
                }
            }
        }
    }

    // Compute m* and find how many moments are realizable
    univariateMomentSet mStar(nMoments_, support_);

    forAll(m0_, celli)
    {
        if (m0_[celli] > 0)
        {
            for (label mi = 0; mi < nMoments_; mi++)
            {
                mStar[mi]
                    = scalar(nFacesOutgoingFlux_[celli] + 1)
                    *moments_[mi][celli] - mPluses[mi][celli];
            }

            nRealizableMomentsStar_[celli] = mStar.nRealizableMoments(false);
        }
        else
        {
            nRealizableMomentsStar_[celli] = nRealizableMoments_[celli];
        }
    }

    // In each cell where the the number of realizable m* is less than the
    // number of realizable m, limitation is attempted
    const cellList& mCells(phi_.mesh().cells());

    forAll(m0_, celli)
    {
        if (nRealizableMomentsStar_[celli] < nRealizableMoments_[celli])
        {
            const cell& mCell(mCells[celli]);

            // Start search for the zetas to limit
            for (label p = 0; p < nRealizableMoments_[celli] - 1; p++)
            {
                scalarList mPlus(nMoments_, 0.0);

                // Check if zeta_p needs limiting by evaluating m* with
                // zeta_k, k > p from constant reconstruction

                // Update mPlus for a face to update m*
                forAll(mCell, fi)
                {
                    const label facei = mCell[fi];

                    if (phi_.mesh().isInternalFace(facei))
                    {
                        if (phi_[facei] > 0)
                        {
                            scalarList zOwn(nZetas_, 0.0);
                            scalarList mOwn(nMoments_, 0.0);

                            for (label zi = 0; zi <= p; zi++)
                            {
                                zOwn[zi] = zetasOwn_[zi][facei];
                            }

                            for (label zi = p + 1; zi < nZetas_; zi++)
                            {
                                zOwn[zi] = zetasUpwindOwn_[zi][facei];
                            }

                            zetaToMoments(zOwn, mOwn, m0Own_[facei]);

                            for (label mi = 0; mi < nMoments_; mi++)
                            {
                                mPlus[mi] += mOwn[mi];
                            }
                        }
                    }
                }

                // Compute m*
                for (label mi = 0; mi < nMoments_; mi++)
                {
                    mStar[mi]
                        = scalar(nFacesOutgoingFlux_[celli] + 1)
                          *moments_[mi][celli] - mPlus[mi];
                }

                nRealizableMomentsStar_[celli]
                    = mStar.nRealizableMoments(false);

                // Check if zeta_p needs limitation
                if (nRealizableMomentsStar_[celli] < nRealizableMoments_[celli])
                {
                    mPlus = 0;

                    // Limit zeta_p
                    forAll(mCell, fi)
                    {
                        const label facei = mCell[fi];

                        if (phi_.mesh().isInternalFace(facei))
                        {
                            if (phi_[facei] > 0)
                            {
                                zetasOwn_[p][facei]
                                    = zetasUpwindOwn_[p][facei]
                                    + 0.5*(zetasCorrOwn_[p][facei]);

                                cellLimiters_[p][celli] = 0.5;

                                scalarList zOwn(nZetas_);
                                scalarList mOwn(nMoments_, 0.0);

                                for (label zi = 0; zi < p; zi++)
                                {
                                    zOwn[zi] = zetasOwn_[zi][facei];
                                }

                                zOwn[p] = zetasOwn_[p][facei];

                                for (label zi = p + 1; zi < nZetas_; zi++)
                                {
                                    zOwn[zi] = zetasUpwindOwn_[zi][facei];
                                }

                                zetaToMoments(zOwn, mOwn, m0Own_[facei]);

                                for (label mi = 0; mi < nMoments_; mi++)
                                {
                                    mPlus[mi] += mOwn[mi];
                                }
                            }
                        }
                    }

                    // Compute m*
                    for (label mi = 0; mi < nMoments_; mi++)
                    {
                        mStar[mi]
                            = scalar(nFacesOutgoingFlux_[celli] + 1)
                              *moments_[mi][celli] - mPlus[mi];
                    }

                    nRealizableMomentsStar_[celli]
                        = mStar.nRealizableMoments(false);

                    if
                    (
                        nRealizableMomentsStar_[celli]
                      < nRealizableMoments_[celli]
                    )
                    {
                        cellLimiters_[p][celli] = 0.0;
                    }
                }
            }
        }
    }

    // Setting limiters on internal faces based on cell limiters
    forAll(phiIf, facei)
    {
        const label own = owner[facei];
        const label nei = neighb[facei];

        if (phi_[facei] > 0)
        {
            for (label zi = 0; zi < nZetas_; zi++)
            {
                limiters_[zi][facei] = cellLimiters_[zi][own];
            }
        }
        else
        {
            for (label zi = 0; zi < nZetas_; zi++)
            {
                limiters_[zi][facei] = cellLimiters_[zi][nei];
            }
        }
    }

    // Setting limiters on boundary faces
    forAll(phiBf, patchi)
    {
        const fvsPatchScalarField& phiPf = phiBf[patchi];

        const labelList& pFaceCells
            = phi_.mesh().boundary()[patchi].faceCells();

        forAll(phiPf, pFacei)
        {
            if (phiPf[pFacei] > 0)
            {
                for (label zi = 0; zi < nZetas_; zi++)
                {
                    limiters_[zi][pFacei]
                        = cellLimiters_[zi][pFaceCells[pFacei]];
                }
            }
        }
    }

    for (label zi = 1; zi < nZetas_; zi++)
    {
        zetasOwn_[zi] = zetasUpwindOwn_[zi] + limiters_[zi]*zetasCorrOwn_[zi];
        zetasNei_[zi] = zetasUpwindNei_[zi] + limiters_[zi]*zetasCorrNei_[zi];
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

    return gMin(internalCo);
}

void Foam::zeta::update()
{
    // Compute zeta fields
    computeZetaFields();

    // Reconstructing zeta_k on cell faces
    interpolateFields();

    // Recompute moments at sides of cell faces
    updateMomentFieldsFromZetas(m0Nei_, zetasNei_, momentsNei_);
    updateMomentFieldsFromZetas(m0Own_, zetasOwn_, momentsOwn_);

    // Apply additional limitation to zeta_k if needed
    limitZetas();

    // Recompute moments at sides of cell faces
    updateMomentFieldsFromZetas(m0Nei_, zetasNei_, momentsNei_);
    updateMomentFieldsFromZetas(m0Own_, zetasOwn_, momentsOwn_);

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
        zetaToMoments(zf, mFace, m0f[facei]);

        for (label mi = 0; mi < nMoments_; mi++)
        {
            mf[mi][facei] = mFace[mi];
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
            zetaToMoments(zf, mFace, m0f.boundaryField()[patchi][facei]);

            for (label mi = 0; mi < nMoments_; mi++)
            {
                mf[mi].boundaryFieldRef()[patchi][facei] = mFace[mi];
            }
        }
    }
}

// ************************************************************************* //
