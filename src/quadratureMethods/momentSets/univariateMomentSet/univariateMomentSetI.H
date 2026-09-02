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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::univariateMomentSet::nZeta() const
{
    return nMoments() - 1;
}


const Foam::scalarList&
Foam::univariateMomentSet::alphaRecurrence() const
{
    return alpha_;
}


Foam::scalarList& Foam::univariateMomentSet::alphaRecurrence()
{
    return alpha_;
}


const Foam::scalarList& Foam::univariateMomentSet::betaRecurrence() const
{
    return beta_;
}


Foam::scalarList& Foam::univariateMomentSet::betaRecurrence()
{
    return beta_;
}


Foam::scalarList& Foam::univariateMomentSet::zetas()
{
    if (!status_.checked)
    {
        checkRealizability(false);
    }

    return zeta_;
}


const Foam::scalarList& Foam::univariateMomentSet::zetas() const
{
    if (!status_.checked)
    {
        FatalErrorInFunction
            << "Zetas are not computed!"
            << exit(FatalError);
    }

    return zeta_;
}


const Foam::supportType& Foam::univariateMomentSet::support() const
{
    return supports_[0];
}


Foam::scalarList& Foam::univariateMomentSet::canonicalMoments()
{
    const supportType& support = supports_[0];

    if (support != supportType::ZeroOne)
    {
        FatalErrorInFunction
            << "Canonical moments are computed only for support [0, 1]." << nl
            << "    The moment set has support "
            << supportTypeToWord(support)
            << "."
            << exit(FatalError);
    }

    if (!status_.checked)
    {
        checkRealizability(false);
    }

    return canonicalMoments_;
}


const Foam::scalarList& Foam::univariateMomentSet::canonicalMoments() const
{
    const supportType& support = supports_[0];

    if (support != supportType::ZeroOne)
    {
        FatalErrorInFunction
            << "Canonical moments are computed only for support [0, 1]." << nl
            << "    The moment set has support " << supportTypeToWord(support)
            << "."
            << exit(FatalError);
    }

    if (!status_.checked)
    {
        FatalErrorInFunction
            << "Canonical moments are not computed!"
            << exit(FatalError);
    }

    return canonicalMoments_;
}


bool Foam::univariateMomentSet::isDegenerate() const
{
    // Checking if moments of order higher than zero are null.
    // The case when m0 = 0 is already managed.
    // Note: the elements are read through the const access operator of the
    // base class, so that the realizability check is not invalidated.
    for (label momenti = 1; momenti < size(); momenti++)
    {
        if (mag(scalarList::operator[](momenti)) > SMALL)
        {
            return false;
        }
    }

    return true;
}


bool Foam::univariateMomentSet::isFullyRealizable
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!status_.checked)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return status_.fullyRealizable;
}


bool Foam::univariateMomentSet::isRealizable
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!status_.checked)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return status_.fullyRealizable || status_.subsetRealizable;
}


bool Foam::univariateMomentSet::isSubsetRealizable
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!status_.checked)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return status_.subsetRealizable;
}


bool Foam::univariateMomentSet::isOnMomentSpaceBoundary
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!status_.checked)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return status_.onMomentSpaceBoundary;
}


Foam::label Foam::univariateMomentSet::nRealizableMoments
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!status_.checked)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return status_.nRealizableMoments;
}


const Foam::scalar&
Foam::univariateMomentSet::moment(const label momentOrder) const
{
    return scalarList::operator[](momentOrder);
}


Foam::scalar& Foam::univariateMomentSet::operator[](const label index)
{
    // Any write to a moment invalidates the quantities derived from the
    // moment vector, and the realizability status
    status_ = realizabilityStatus();

    return scalarList::operator[](index);
}


const Foam::scalar&
Foam::univariateMomentSet::operator[](const label index) const
{
    return scalarList::operator[](index);
}


const Foam::scalar&
Foam::univariateMomentSet::operator()(const label momentOrder) const
{
    return scalarList::operator[](momentOrder);
}


Foam::scalar& Foam::univariateMomentSet::operator()(const label momentOrder)
{
    return operator[](momentOrder);
}


// ************************************************************************* //
