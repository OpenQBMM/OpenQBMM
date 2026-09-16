/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2026 Alberto Passalacqua
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

#include "zetaUnivariateAdvection.H"
#include "upwind.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace univariateAdvection
{
    defineTypeNameAndDebug(zeta, 0);

    addToRunTimeSelectionTable
    (
        univariateMomentAdvection,
        zeta,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateAdvection::zeta::zeta
(
    const dictionary& dict,
    const scalarQuadratureApproximation& quadrature,
    const surfaceScalarField& phi,
    const supportType& support
)
:
    univariateMomentAdvection(dict, quadrature, phi, support),
    m0_(moments_(0)),
    m0Own_
    (
        IOobject::groupName("m0OwnZeta", name_),
        fvc::interpolate(m0_, own_, "reconstruct(m0)")
    ),
    m0Nei_
    (
        IOobject::groupName("m0NeiZeta", name_),
        fvc::interpolate(m0_, nei_, "reconstruct(m0)")
    ),
    nAuxiliaryFields_(nMoments_ - 1),
    auxiliaryFields_(nAuxiliaryFields_),
    auxiliaryFieldsNei_(nAuxiliaryFields_),
    auxiliaryFieldsOwn_(nAuxiliaryFields_),
    auxiliaryFieldsUpwindNei_(nAuxiliaryFields_),
    auxiliaryFieldsUpwindOwn_(nAuxiliaryFields_),
    auxiliaryFieldsCorrNei_(nAuxiliaryFields_),
    auxiliaryFieldsCorrOwn_(nAuxiliaryFields_),
    momentsNei_(nMoments_),
    momentsOwn_(nMoments_),
    nFacesOutgoingFlux_(m0_.size(), 0),
    nRealizableMoments_(m0_.size(), 0),
    nRealizableMomentsStar_(m0_.size(), 0),
    limiters_(nAuxiliaryFields_),
    cellLimiters_(nAuxiliaryFields_),
    smallM0_(dict.lookupOrDefault<scalar>("smallM0", SMALL)),
    smallZeta_(dict.lookupOrDefault<scalar>("smallZeta", SMALL)),
    smallAuxiliaryQuantity_
    (
        dict.lookupOrDefault<scalar>("smallAuxiliaryQuantity", 1.0e-7)
    )
{
    if
    (
        quadrature.momentOrders()[0].size() > 1
     || (support_ != supportType::RPlus && support_ != supportType::ZeroOne)
    )
    {
        FatalErrorInFunction
            << "Zeta advection scheme can only be used for" << nl
            << "univariate distributions with support over R+ or [0, 1]"
            << abort(FatalError);
    }

    if (support_ == supportType::RPlus)
    {
        Info << endl << "Using zeta scheme with R+ support.\n" << endl;
    }
    else
    {
        Info << endl << "Using zeta scheme with [0, 1] support.\n" << endl;
    }

    // Populating zeta_k fields and interpolated zeta_k fields
    forAll(auxiliaryFields_, auxiliaryFieldi)
    {
        auxiliaryFields_.set
        (
            auxiliaryFieldi,
            new volScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryField", {auxiliaryFieldi}),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        auxiliaryFieldsNei_.set
        (
            auxiliaryFieldi,
            new surfaceScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryFieldNei", {auxiliaryFieldi}),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        auxiliaryFieldsOwn_.set
        (
            auxiliaryFieldi,
            new surfaceScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryFieldOwn", {auxiliaryFieldi}),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        auxiliaryFieldsUpwindNei_.set
        (
            auxiliaryFieldi,
            new surfaceScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryFieldUpwindNei", {auxiliaryFieldi}),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        auxiliaryFieldsUpwindOwn_.set
        (
            auxiliaryFieldi,
            new surfaceScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryFieldUpwindOwn", {auxiliaryFieldi}),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        auxiliaryFieldsCorrNei_.set
        (
            auxiliaryFieldi,
            new surfaceScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryFieldCorrNei", {auxiliaryFieldi}),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        auxiliaryFieldsCorrOwn_.set
        (
            auxiliaryFieldi,
            new surfaceScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryFieldCorrOwn", {auxiliaryFieldi}),
                    phi.mesh().time().timeName(),
                    phi.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                phi.mesh(),
                dimensionedScalar("zero", dimless, Zero)
            )
        );

        limiters_.set
        (
            auxiliaryFieldi,
            new surfaceScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryFieldLimiter", {auxiliaryFieldi}),
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
            auxiliaryFieldi,
            new volScalarField
            (
                IOobject
                (
                    fieldName("auxiliaryFieldCellLimiter", {auxiliaryFieldi}),
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
                fieldName("momentNeiZeta", {momenti}),
                fvc::interpolate(moments_(momenti))
            )
        );

        momentsOwn_.set
        (
            momenti,
            new surfaceScalarField
            (
                fieldName("momentOwnZeta", {momenti}),
                fvc::interpolate(moments_(momenti))
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateAdvection::zeta::~zeta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::univariateAdvection::zeta::interpolateFields()
{
    // The limiter is part of the scheme rather than a choice left to the
    // case: the second-order reconstruction is realizable because it is
    // limited this way, so it is not read from the dictionary
    IStringStream m0OwnLimiter("Minmod");
    IStringStream auxiliaryFieldsOwnLimiter("Minmod");

    tmp<surfaceInterpolationScheme<scalar>> m0OwnScheme
    (
        fvc::scheme<scalar>(own_, m0OwnLimiter)
    );

    tmp<surfaceInterpolationScheme<scalar>> auxiliaryFieldsOwnScheme
    (
        fvc::scheme<scalar>(own_, auxiliaryFieldsOwnLimiter)
    );

    IStringStream m0NeiLimiter("Minmod");
    IStringStream auxiliaryFieldsNeiLimiter("Minmod");

    tmp<surfaceInterpolationScheme<scalar>> m0NeiScheme
    (
        fvc::scheme<scalar>(nei_, m0NeiLimiter)
    );

    tmp<surfaceInterpolationScheme<scalar>> auxiliaryFieldsNeiScheme
    (
        fvc::scheme<scalar>(nei_, auxiliaryFieldsNeiLimiter)
    );

    m0Own_ = m0OwnScheme().interpolate(moments_(0));
    m0Nei_ = m0NeiScheme().interpolate(moments_(0));

    forAll(auxiliaryFields_, fieldi)
    {
        auxiliaryFieldsNei_[fieldi] =
            auxiliaryFieldsNeiScheme().interpolate(auxiliaryFields_[fieldi]);

        auxiliaryFieldsOwn_[fieldi] =
            auxiliaryFieldsOwnScheme().interpolate(auxiliaryFields_[fieldi]);

        // The upwind value of each side is the value of the cell on that
        // side. This used to be the flux of the upwind scheme, which is
        // the interpolated value times the flux the scheme is built with:
        // one on the owner side, but minus one on the neighbour side, so
        // the neighbour side carried the opposite of its cell value, and
        // any limiter below one reconstructed a negative auxiliary
        // quantity there.
        auxiliaryFieldsUpwindNei_[fieldi] =
            upwind<scalar>
            (
                auxiliaryFields_[fieldi].mesh(), nei_
            ).interpolate(auxiliaryFields_[fieldi]);

        auxiliaryFieldsUpwindOwn_[fieldi] =
            upwind<scalar>
            (
                auxiliaryFields_[fieldi].mesh(), own_
            ).interpolate(auxiliaryFields_[fieldi]);

        auxiliaryFieldsCorrNei_[fieldi] =
            auxiliaryFieldsNei_[fieldi] - auxiliaryFieldsUpwindNei_[fieldi];

        auxiliaryFieldsCorrOwn_[fieldi] =
            auxiliaryFieldsOwn_[fieldi] - auxiliaryFieldsUpwindOwn_[fieldi];
    }
}

void Foam::univariateAdvection::zeta::zetaToMoments
(
    const scalarList& zetaf,
    scalarList& mf,
    const label nMoments,
    scalar m0
)
{
    // A set of nMoments moments is described by its zero-order moment and
    // the nMoments - 1 values of the zeta chain, so all of them are needed
    // to recover it. Both recursions below therefore run to nMoments: the
    // index of zeta they reach is j - i and i - 1 respectively, which stays
    // within the chain. Stopping them one short, as was done previously,
    // leaves the last value of zeta unread and the moment of highest order
    // wrong - two nodes of equal weight at 1 and 2 returned a third moment
    // of 25/6 rather than 9/2.
    scalarSquareMatrix S(nMoments, 0.0);

    for (label i = 0; i < nMoments; i++)
    {
        S[0][i] = 1.0;
    }

    for (label i = 1; i < nMoments; i++)
    {
        for (label j = i; j < nMoments; j++)
        {
            S[i][j] = S[i][j - 1] + zetaf[j - i]*S[i - 1][j];
        }
    }

    // prod[i] is the product of the first i values of the zeta chain
    scalarList prod(nMoments, 1.0);

    prod[1] = zetaf[0];

    for (label i = 2; i < nMoments; i++)
    {
        prod[i] = prod[i - 1]*zetaf[i - 1];
    }

    // Resetting moments to zero
    mf = 0.0;

    // Computing moments
    mf[0] = 1.0;
    mf[1] = zetaf[0];

    for (label i = 2; i < nMoments; i++)
    {
        for (label j = 0; j <= i/2; j++)
        {
            mf[i] += prod[i - 2*j]*sqr(S[j][i - j]);
        }
    }

    if (m0 != 1.0)
    {
        for (label mi = 0; mi < nMoments; mi++)
        {
            mf[mi] *= m0;
        }
    }
}

void Foam::univariateAdvection::zeta::canonicalMomentsToMoments
(
    const scalarList& canonicalMomentsf,
    scalarList& mf,
    const label nMoments,
    scalar m0
)
{
    const label nAuxiliaryFields = nMoments - 1;

    scalarList zetas(nAuxiliaryFields);
    zetas[0] = canonicalMomentsf[0];

    for (label i = 1; i < nAuxiliaryFields; i++)
    {
        zetas[i] = canonicalMomentsf[i]*(1.0 - canonicalMomentsf[i - 1]);
    }

    zetaToMoments(zetas, mf, nMoments, m0);
}

void Foam::univariateAdvection::zeta::auxiliaryQuantitiesToMoments
(
    const scalarList& auxiliaryQuantityf,
    scalarList& mf,
    scalar m0
)
{
    if (support_ == supportType::RPlus)
    {
        zetaToMoments(auxiliaryQuantityf, mf, nMoments_, m0);
    }
    else // Support is [0, 1]
    {
        canonicalMomentsToMoments(auxiliaryQuantityf, mf, nMoments_, m0);
    }
}

bool Foam::univariateAdvection::zeta::outgoingFace
(
    const label celli,
    const label facei,
    label& patchi,
    label& pFacei,
    bool& ownSide
) const
{
    const fvMesh& mesh = phi_.mesh();

    if (mesh.isInternalFace(facei))
    {
        patchi = -1;
        pFacei = facei;

        // The flux leaves through the owner side of a face the cell owns
        // and through the neighbour side of a face it is the neighbour of
        if (phi_[facei] > 0)
        {
            ownSide = true;

            return mesh.owner()[facei] == celli;
        }
        else if (phi_[facei] < 0)
        {
            ownSide = false;

            return mesh.neighbour()[facei] == celli;
        }

        return false;
    }

    patchi = mesh.boundaryMesh().whichPatch(facei);

    if (patchi < 0)
    {
        return false;
    }

    pFacei = facei - mesh.boundaryMesh()[patchi].start();
    ownSide = true;

    const surfaceScalarField::Boundary& phiBf = phi_.boundaryField();

    // Patches without a finite volume representation, such as empty and
    // wedge patches, carry no flux
    if (pFacei >= phiBf[patchi].size())
    {
        return false;
    }

    return phiBf[patchi][pFacei] > 0;
}


void Foam::univariateAdvection::zeta::addFaceMomentsToMPlus
(
    const label p,
    const label patchi,
    const label facei,
    const bool ownSide,
    scalarList& mPlus
)
{
    const PtrList<surfaceScalarField>& limited =
        ownSide ? auxiliaryFieldsOwn_ : auxiliaryFieldsNei_;

    const PtrList<surfaceScalarField>& upwind =
        ownSide ? auxiliaryFieldsUpwindOwn_ : auxiliaryFieldsUpwindNei_;

    const surfaceScalarField& m0f = ownSide ? m0Own_ : m0Nei_;

    scalarList auxiliaryQuantity(nAuxiliaryFields_, Zero);
    scalarList mFace(nMoments_, Zero);
    scalar m0Face = 0.0;

    if (patchi < 0)
    {
        for (label i = 0; i <= p; i++)
        {
            auxiliaryQuantity[i] = limited[i][facei];
        }

        for (label i = p + 1; i < nAuxiliaryFields_; i++)
        {
            auxiliaryQuantity[i] = upwind[i][facei];
        }

        m0Face = m0f[facei];
    }
    else
    {
        for (label i = 0; i <= p; i++)
        {
            auxiliaryQuantity[i] = limited[i].boundaryField()[patchi][facei];
        }

        for (label i = p + 1; i < nAuxiliaryFields_; i++)
        {
            auxiliaryQuantity[i] = upwind[i].boundaryField()[patchi][facei];
        }

        m0Face = m0f.boundaryField()[patchi][facei];
    }

    auxiliaryQuantitiesToMoments(auxiliaryQuantity, mFace, m0Face);

    for (label mi = 0; mi < nMoments_; mi++)
    {
        mPlus[mi] += mFace[mi];
    }
}


void Foam::univariateAdvection::zeta::computeAuxiliaryFields()
{
    // The moment set is allocated once and reused over cells and faces, so
    // that its lists and index map are not rebuilt for every cell
    univariateMomentSet m(nMoments_, support_, smallM0_, smallZeta_);

    // Cell-center values
    forAll(m0_, celli)
    {
        if (m0_[celli] < smallM0_)
        {
            // A cell holding no distribution has no auxiliary quantities.
            // They used to be left at the values of the previous timestep,
            // along with the count of realizable moments, and the limiter
            // then worked against both.
            for (label i = 0; i < nAuxiliaryFields_; i++)
            {
                auxiliaryFields_[i][celli] = 0.0;
            }

            nRealizableMoments_[celli] = 0;
        }
        else
        {
            for (label mi = 0; mi < nMoments_; mi++)
            {
                m[mi] = moments_(mi)[celli];
            }

            nRealizableMoments_[celli] = m.nRealizableMoments();

            // AP: Recover zeta_k if support is R+, otherwise obtain
            //     canonical moments. Support over R is excluded in the
            //     constructor so it cannot be encountered.
            scalarList& auxiliaryQuantities
            (
                support_ == supportType::RPlus ? m.zetas() : m.canonicalMoments()
            );

            for (label i = 0; i < nAuxiliaryFields_; i++)
            {
                // Both arms of the test this replaces assigned the same
                // thing, so what it does is clip to zero from below
                auxiliaryFields_[i][celli] =
                    auxiliaryQuantities[i] > smallAuxiliaryQuantity_
                  ? auxiliaryQuantities[i]
                  : 0.0;
            }
        }
    }

    // Boundary conditions
    const volScalarField::Boundary& bf = auxiliaryFields_[0].boundaryField();

    forAll(bf, patchi)
    {
        const fvPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            if (m0_.boundaryField()[patchi][facei] < smallM0_)
            {
                // A face holding no distribution has no auxiliary
                // quantities, like a cell; this tested a bare SMALL where
                // the cells test smallM0, and left the face as it was
                for (label i = 0; i < nAuxiliaryFields_; i++)
                {
                    auxiliaryFields_[i].boundaryFieldRef()[patchi][facei] = 0.0;
                }
            }
            else
            {
                for (label mi = 0; mi < nMoments_; mi++)
                {
                    m[mi] = moments_(mi).boundaryField()[patchi][facei];
                }

                // AP: Recover zeta_k if support is R+, otherwise obtain
                //     canonical moments. Support over R is excluded in the
                //     constructor so it cannot be encountered.
                scalarList& auxiliaryQuantities
                (
                    support_ == supportType::RPlus ? m.zetas() : m.canonicalMoments()
                );

                for (label i = 0; i < nAuxiliaryFields_; i++)
                {
                    volScalarField& auxiliaryFieldi = auxiliaryFields_[i];

                    volScalarField::Boundary& auxiliaryFieldiBf =
                        auxiliaryFieldi.boundaryFieldRef();

                    auxiliaryFieldiBf[patchi][facei] = auxiliaryQuantities[i];
                }
            }
        }
    }

    forAll(auxiliaryFields_, fieldi)
    {
        auxiliaryFields_[fieldi].correctBoundaryConditions();
    }
}

void Foam::univariateAdvection::zeta::countFacesWithOutgoingFlux() const
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

void Foam::univariateAdvection::zeta::limitAuxiliaryFields()
{
    const fvMesh& mesh = phi_.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const scalarField& phiIf = phi_;
    const surfaceScalarField::Boundary& phiBf = phi_.boundaryField();
    const label nInternalFaces = mesh.nInternalFaces();

    countFacesWithOutgoingFlux();

    forAll(cellLimiters_, li)
    {
        forAll(cellLimiters_[0], celli)
        {
            cellLimiters_[li][celli] = 1.0;
        }

        // The face limiters are reset too: a boundary face the flux does
        // not leave through is not set below, and kept the limiter of a
        // step at which the flux did
        limiters_[li] = dimensionedScalar(dimless, 1.0);
    }

    // First check on m* to identify cells in need of additional limitation
    scalarRectangularMatrix mPluses(nMoments_, m0_.size(), 0.0);

    // Find m+ (moments reconstructed on cell faces with outgoing flux)
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        if (phi_[facei] > 0.0)
        {
            for (label mi = 0; mi < nMoments_; mi++)
            {
                mPluses[mi][own] += momentsOwn_[mi][facei];
            }
        }
        else if (phi_[facei] < 0.0)
        {
            // A face without flux is not counted by
            // countFacesWithOutgoingFlux, so it does not add to m+ either;
            // it used to, biasing m* of the neighbour towards unrealizable
            // on every face the flow is tangent to
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
                    mPluses[mi][pFaceCells[pFacei]] +=
                        momentsOwn_[mi].boundaryField()[patchi][pFacei];
                }
            }
        }
    }

    // Compute m* and find how many moments are realizable
    univariateMomentSet mStar(nMoments_, support_, smallM0_, smallZeta_);

    forAll(m0_, celli)
    {
        if (m0_[celli] >= smallM0_)
        {
            for (label mi = 0; mi < nMoments_; mi++)
            {
                mStar[mi]
                    = scalar(nFacesOutgoingFlux_[celli] + 1)
                        *moments_(mi)[celli] - mPluses[mi][celli];
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

            // Start search for the auxiliary quantities to limit
            for (label p = 0; p < nRealizableMoments_[celli] - 1; p++)
            {
                scalarList mPlus(nMoments_, Zero);

                // Check if the auxiliary quantity with index p needs limiting
                // by evaluating m* with auxiliaryQuantity_k, k > p from
                // constant reconstruction. m+ is the sum over the faces the
                // flux leaves the cell through, each reconstructed on the
                // side of it the cell is on: the search used to take every
                // face of the cell with a positive flux, which is outgoing
                // only for the cell that owns it, and to leave out the
                // faces the cell is the neighbour of, so m* was wrong in
                // any cell with a face of either kind. Boundary faces count
                // because countFacesWithOutgoingFlux includes them.
                forAll(mCell, fi)
                {
                    label patchi = -1;
                    label pFacei = -1;
                    bool ownSide = true;

                    if (outgoingFace(celli, mCell[fi], patchi, pFacei, ownSide))
                    {
                        addFaceMomentsToMPlus(p, patchi, pFacei, ownSide, mPlus);
                    }
                }

                // Compute m*
                for (label mi = 0; mi < nMoments_; mi++)
                {
                    mStar[mi]
                        = scalar(nFacesOutgoingFlux_[celli] + 1)
                          *moments_(mi)[celli] - mPlus[mi];
                }

                nRealizableMomentsStar_[celli]
                    = mStar.nRealizableMoments(false);

                // Check if the auxiliary quantity of index p needs limitation
                if (nRealizableMomentsStar_[celli] < nRealizableMoments_[celli])
                {
                    mPlus = 0;

                    // Limit the auxiliary quantity of index p to half of
                    // its correction on every face the flux leaves through,
                    // on the side the cell reconstructs
                    forAll(mCell, fi)
                    {
                        label patchi = -1;
                        label pFacei = -1;
                        bool ownSide = true;

                        if
                        (
                            !outgoingFace
                            (
                                celli, mCell[fi], patchi, pFacei, ownSide
                            )
                        )
                        {
                            continue;
                        }

                        surfaceScalarField& limited =
                            ownSide
                          ? auxiliaryFieldsOwn_[p]
                          : auxiliaryFieldsNei_[p];

                        const surfaceScalarField& upwind =
                            ownSide
                          ? auxiliaryFieldsUpwindOwn_[p]
                          : auxiliaryFieldsUpwindNei_[p];

                        const surfaceScalarField& correction =
                            ownSide
                          ? auxiliaryFieldsCorrOwn_[p]
                          : auxiliaryFieldsCorrNei_[p];

                        if (patchi < 0)
                        {
                            limited[pFacei] =
                                upwind[pFacei] + 0.5*correction[pFacei];
                        }
                        else
                        {
                            limited.boundaryFieldRef()[patchi][pFacei] =
                                upwind.boundaryField()[patchi][pFacei]
                              + 0.5*correction.boundaryField()[patchi][pFacei];
                        }

                        cellLimiters_[p][celli] = 0.5;

                        addFaceMomentsToMPlus(p, patchi, pFacei, ownSide, mPlus);
                    }

                    // Compute m*
                    for (label mi = 0; mi < nMoments_; mi++)
                    {
                        mStar[mi] =
                            scalar(nFacesOutgoingFlux_[celli] + 1)
                           *moments_(mi)[celli] - mPlus[mi];
                    }

                    nRealizableMomentsStar_[celli] =
                        mStar.nRealizableMoments(false);

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
        const label nei = neighbour[facei];

        if (phi_[facei] > 0)
        {
            for (label i = 0; i < nAuxiliaryFields_; i++)
            {
                limiters_[i][facei] = cellLimiters_[i][own];
            }
        }
        else
        {
            for (label i = 0; i < nAuxiliaryFields_; i++)
            {
                limiters_[i][facei] = cellLimiters_[i][nei];
            }
        }
    }

    // Setting limiters on boundary faces. A face the flux leaves through
    // takes the limiter of its cell. A face of a coupled patch the flux
    // enters through takes the limiter of the cell across it, so that the
    // two sides of one face reconstruct the same moments: left at one on
    // the entering side while the leaving side was limited, the flux out
    // of a cell through a cyclic or processor boundary was not the flux
    // into the cell across it, and the moments were not conserved.
    forAll(cellLimiters_, i)
    {
        cellLimiters_[i].correctBoundaryConditions();
    }

    forAll(phiBf, patchi)
    {
        const fvsPatchScalarField& phiPf = phiBf[patchi];
        const fvPatch& patch = mesh.boundary()[patchi];
        const labelList& pFaceCells = patch.faceCells();

        for (label i = 0; i < nAuxiliaryFields_; i++)
        {
            scalarField& limiterPf = limiters_[i].boundaryFieldRef()[patchi];

            tmp<scalarField> tacross;

            if (patch.coupled())
            {
                tacross =
                    cellLimiters_[i].boundaryField()[patchi]
                   .patchNeighbourField();
            }

            forAll(phiPf, pFacei)
            {
                if (phiPf[pFacei] > 0)
                {
                    limiterPf[pFacei] = cellLimiters_[i][pFaceCells[pFacei]];
                }
                else if (patch.coupled() && phiPf[pFacei] < 0)
                {
                    limiterPf[pFacei] = tacross()[pFacei];
                }
            }
        }
    }

    // The first auxiliary quantity is m1/m0, on which the first moment
    // depends linearly, so it is limited like all the others. The search
    // above writes a trial value of 0.5 into auxiliaryFieldsOwn_ as it goes;
    // reconstructing from index 0 is what replaces that trial value by the
    // limiter the search settled on, and is the only thing that limits
    // auxiliaryFieldsNei_ at all.
    for (label i = 0; i < nAuxiliaryFields_; i++)
    {
        auxiliaryFieldsOwn_[i] =
            auxiliaryFieldsUpwindOwn_[i]
          + limiters_[i]*auxiliaryFieldsCorrOwn_[i];

        auxiliaryFieldsNei_[i] =
            auxiliaryFieldsUpwindNei_[i]
          + limiters_[i]*auxiliaryFieldsCorrNei_[i];
    }
}

Foam::scalar Foam::univariateAdvection::zeta::realizableCo() const
{
    // The realizability condition of the scheme is written in terms of the
    // number of faces of a cell carrying an outgoing flux, so the same count
    // the limiter works with is used here. Counting only internal faces, as
    // was done previously, returns a Courant limit a boundary cell does not
    // actually satisfy.
    countFacesWithOutgoingFlux();

    scalarField co(m0_.size(), Zero);

    forAll(co, celli)
    {
        co[celli] = 1.0/scalar(nFacesOutgoingFlux_[celli] + 1);
    }

    return gMin(co);
}

void Foam::univariateAdvection::zeta::update()
{
    if (m0_.size() != nFacesOutgoingFlux_.size())
    {
        nFacesOutgoingFlux_.resize(m0_.size(), 0);
        nRealizableMoments_.resize(m0_.size(), 0);
        nRealizableMomentsStar_.resize(m0_.size(), 0);
    }

    // Compute zeta fields
    computeAuxiliaryFields();

    // Reconstructing auxiliary fields on cell faces
    interpolateFields();

    // Recompute moments at sides of cell faces
    updateMomentFieldsFromAuxiliaryQuantities
    (
        m0Nei_, auxiliaryFieldsNei_, momentsNei_
    );

    updateMomentFieldsFromAuxiliaryQuantities
    (
        m0Own_, auxiliaryFieldsOwn_, momentsOwn_
    );

    // Apply additional limitation to auxiliary quantities if needed
    limitAuxiliaryFields();

    // Recompute moments at sides of cell faces
    updateMomentFieldsFromAuxiliaryQuantities
    (
        m0Nei_, auxiliaryFieldsNei_, momentsNei_
    );

    updateMomentFieldsFromAuxiliaryQuantities
    (
        m0Own_, auxiliaryFieldsOwn_, momentsOwn_
    );

    // Calculate moment advection term
    dimensionedScalar zeroPhi("zero", phi_.dimensions(), Zero);

    forAll(divMoments_, divi)
    {
        divMoments_(divi) =
            fvc::surfaceIntegrate
            (
                momentsNei_[divi]*min(phi_, zeroPhi)
              + momentsOwn_[divi]*max(phi_, zeroPhi)
            );
    }
}

void Foam::univariateAdvection::zeta::updateMomentFieldsFromAuxiliaryQuantities
(
    const surfaceScalarField& m0f,
    const PtrList<surfaceScalarField>& auxiliaryFieldsf,
    PtrList<surfaceScalarField>& mf
)
{
    forAll(auxiliaryFieldsf[0], facei)
    {
        scalarList auxiliaryQuantitiesf(nAuxiliaryFields_);

        for (label i = 0; i < nAuxiliaryFields_; i++)
        {
            auxiliaryQuantitiesf[i] = auxiliaryFieldsf[i][facei];
        }

        scalarList mFace(nMoments_, Zero);
        auxiliaryQuantitiesToMoments(auxiliaryQuantitiesf, mFace, m0f[facei]);

        for (label mi = 0; mi < nMoments_; mi++)
        {
            mf[mi][facei] = mFace[mi];
        }
    }

    // Boundary conditions
    const surfaceScalarField::Boundary& bf =
        auxiliaryFieldsf[0].boundaryField();

    forAll(bf, patchi)
    {
        const fvsPatchScalarField& m0Patch = bf[patchi];

        forAll(m0Patch, facei)
        {
            scalarList auxiliaryQuantitiesf(nAuxiliaryFields_);

            for (label i = 0; i < nAuxiliaryFields_; i++)
            {
                auxiliaryQuantitiesf[i] =
                    auxiliaryFieldsf[i].boundaryField()[patchi][facei];
            }

            scalarList mFace(nMoments_, Zero);

            auxiliaryQuantitiesToMoments
            (
                auxiliaryQuantitiesf, mFace, m0f.boundaryField()[patchi][facei]
            );

            for (label mi = 0; mi < nMoments_; mi++)
            {
                mf[mi].boundaryFieldRef()[patchi][facei] = mFace[mi];
            }
        }
    }
}

// ************************************************************************* //
