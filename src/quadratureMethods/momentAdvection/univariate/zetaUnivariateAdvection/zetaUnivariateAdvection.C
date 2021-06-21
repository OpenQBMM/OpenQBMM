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
    Copyright (C) 2019-2021 Alberto Passalacqua
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
    const word& support
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
    phi_(phi)
{
    if 
    (
        quadrature.momentOrders()[0].size() > 1
     || (support_ != "RPlus" && support_ != "01")
    )
    {
        FatalErrorInFunction
            << "Zeta advection scheme can only be used for" << nl
            << "univariate distributions with support over R+ or [0, 1]"
            << abort(FatalError);
    }

    if (support_ == "RPlus")
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
                    IOobject::AUTO_WRITE
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

        auxiliaryFieldsUpwindNei_[fieldi] =
            upwind<scalar>
            (
                auxiliaryFields_[fieldi].mesh(), nei_
            ).flux(auxiliaryFields_[fieldi]);

        auxiliaryFieldsUpwindOwn_[fieldi] =
            upwind<scalar>
            (
                auxiliaryFields_[fieldi].mesh(), own_
            ).flux(auxiliaryFields_[fieldi]);

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
    scalar m0
)
{
    scalarSquareMatrix S(nMoments_, 0.0);

    for (label i = 0; i < nAuxiliaryFields_; i++)
    {
        S[0][i] = 1.0;
    }

    for (label i = 1; i < nAuxiliaryFields_; i++)
    {
        for (label j = i; j < nAuxiliaryFields_; j++)
        {
            S[i][j] = S[i][j - 1] + zetaf[j - i]*S[i - 1][j];
        }
    }

    scalarList prod(nMoments_, 1.0);

    prod[1] = zetaf[0];

    for (label i = 2; i < nAuxiliaryFields_; i++)
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

void Foam::univariateAdvection::zeta::canonicalMomentsToMoments
(
    const scalarList& canonicalMomentsf,
    scalarList& mf,
    scalar m0
)
{
    scalarList zetas(nAuxiliaryFields_);
    zetas[0] = canonicalMomentsf[0];
    
    for (label i = 1; i < nAuxiliaryFields_; i++)
    {
        zetas[i] = canonicalMomentsf[i]*(1.0 - canonicalMomentsf[i - 1]);
    }

    zetaToMoments(zetas, mf, m0);
}

void Foam::univariateAdvection::zeta::auxiliaryQuantitiesToMoments
(
    const scalarList& auxiliaryQuantityf,
    scalarList& mf,
    scalar m0
)
{
    if (support_ == "RPlus")
    {
        zetaToMoments(auxiliaryQuantityf, mf, m0);
    }
    else // Support is [0, 1]
    {
        canonicalMomentsToMoments(auxiliaryQuantityf, mf, m0);
    }
}

void Foam::univariateAdvection::zeta::computeAuxiliaryFields()
{
    // Cell-center values
    forAll(m0_, celli)
    {
        if (m0_[celli] >= SMALL)
        {
            univariateMomentSet m(nMoments_, support_);

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
                support_ == "RPlus" ? m.zetas() : m.canonicalMoments()
            );

            for (label i = 0; i < nAuxiliaryFields_; i++)
            {
                auxiliaryFields_[i][celli] = auxiliaryQuantities[i];

                if (auxiliaryFields_[i][celli] > 1.0e-7)
                {
                    auxiliaryFields_[i][celli] = auxiliaryQuantities[i];
                }
                else
                {
                    auxiliaryFields_[i][celli] = 0.0;
                }
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
            if (m0_.boundaryField()[patchi][facei] >= SMALL)
            {
                univariateMomentSet m(nMoments_, support_);

                for (label mi = 0; mi < nMoments_; mi++)
                {
                    m[mi] = moments_(mi).boundaryField()[patchi][facei];
                }

                // AP: Recover zeta_k if support is R+, otherwise obtain
                //     canonical moments. Support over R is excluded in the
                //     constructor so it cannot be encountered.
                scalarList& auxiliaryQuantities
                (
                    support_ == "RPlus" ? m.zetas() : m.canonicalMoments()
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
}

void Foam::univariateAdvection::zeta::countFacesWithOutgoingFlux()
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
    const labelUList& owner = phi_.mesh().owner();
    const labelUList& neighbour = phi_.mesh().neighbour();
    const scalarField& phiIf = phi_;
    const surfaceScalarField::Boundary& phiBf = phi_.boundaryField();
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
                // constant reconstruction

                // Update mPlus for a face to update m*
                forAll(mCell, fi)
                {
                    const label facei = mCell[fi];

                    if (phi_.mesh().isInternalFace(facei))
                    {
                        if (phi_[facei] > 0)
                        {
                            scalarList auxiliaryQuantityOwn
                            (
                                nAuxiliaryFields_, 
                                Zero
                            );

                            scalarList mOwn(nMoments_, Zero);

                            for (label i = 0; i <= p; i++)
                            {
                                auxiliaryQuantityOwn[i] = 
                                    auxiliaryFieldsOwn_[i][facei];
                            }

                            for (label i = p + 1; i < nAuxiliaryFields_; i++)
                            {
                                auxiliaryQuantityOwn[i] =
                                    auxiliaryFieldsUpwindOwn_[i][facei];
                            }

                            auxiliaryQuantitiesToMoments(
                                auxiliaryQuantityOwn, mOwn, m0Own_[facei]);

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
                          *moments_(mi)[celli] - mPlus[mi];
                }

                nRealizableMomentsStar_[celli]
                    = mStar.nRealizableMoments(false);

                // Check if the auxiliary quantity of index p needs limitation
                if (nRealizableMomentsStar_[celli] < nRealizableMoments_[celli])
                {
                    mPlus = 0;

                    // Limit auxiliary quantities
                    forAll(mCell, fi)
                    {
                        const label facei = mCell[fi];

                        if (phi_.mesh().isInternalFace(facei))
                        {
                            if (phi_[facei] > 0)
                            {
                                auxiliaryFieldsOwn_[p][facei] = 
                                    auxiliaryFieldsUpwindOwn_[p][facei]
                                  + 0.5*(auxiliaryFieldsCorrOwn_[p][facei]);

                                cellLimiters_[p][celli] = 0.5;

                                scalarList auxiliaryQuantityOwn
                                (
                                    nAuxiliaryFields_
                                );

                                scalarList mOwn(nMoments_, Zero);

                                for (label i = 0; i < p; i++)
                                {
                                    auxiliaryQuantityOwn[i] = 
                                        auxiliaryFieldsOwn_[i][facei];
                                }

                                auxiliaryQuantityOwn[p] = 
                                    auxiliaryFieldsOwn_[p][facei];

                                for 
                                (
                                    label i = p + 1; 
                                    i < nAuxiliaryFields_; 
                                    i++)
                                {
                                    auxiliaryQuantityOwn[i] = 
                                        auxiliaryFieldsUpwindOwn_[i][facei];
                                }

                                auxiliaryQuantitiesToMoments
                                (
                                    auxiliaryQuantityOwn, mOwn, m0Own_[facei]
                                );

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
                for (label i = 0; i < nAuxiliaryFields_; i++)
                {
                    limiters_[i][pFacei] = cellLimiters_[i][pFaceCells[pFacei]];
                }
            }
        }
    }

    for (label i = 1; i < nAuxiliaryFields_; i++)
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
    const fvMesh& mesh(phi_.mesh());
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    scalarField internalCo(m0_.size(), Zero);

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

void Foam::univariateAdvection::zeta::update()
{
    if (m0_.size() != nFacesOutgoingFlux_.size())
    {
        nFacesOutgoingFlux_.resize(m0_.size());
        nRealizableMoments_.resize(m0_.size());
        nRealizableMomentsStar_.resize(m0_.size());
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