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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

Foam::label Foam::univariateMomentSet::negativeZeta()
{
    if (!realizabilityChecked_)
    {
        checkRealizability(false);
    }

    return negativeZeta_;
}

Foam::scalarList& Foam::univariateMomentSet::zetas()
{
    if (!realizabilityChecked_)
    {
        checkRealizability(false);
    }

    return zeta_;
}

Foam::scalarList& Foam::univariateMomentSet::canonicalMoments()
{
    if (support_ != "01")
    {
        FatalErrorInFunction
            << "Canonical moments are computed only for support [0, 1]." << nl
            << "    The moment set has support " << support_ << "."
            << abort(FatalError);
    }

    if (!realizabilityChecked_)
    {
        checkRealizability(false);
    }

    return canonicalMoments_;
}

bool Foam::univariateMomentSet::isDegenerate()
{
    // Checking if moments of order higher than zero are null
    // The case when m0 = 0 is already managed
    for (Foam::label momentI = 1; momentI < (*this).size(); momentI++)
    {
        if (mag((*this)[momentI]) > SMALL)
        {
            degenerate_ = false;

            return degenerate_;
        }
    }

    degenerate_ = true;

    return degenerate_;
}

bool Foam::univariateMomentSet::isFullyRealizable
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!realizabilityChecked_)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return fullyRealizable_;
}

bool Foam::univariateMomentSet::isRealizable
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!realizabilityChecked_)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    bool realizable = fullyRealizable_ || subsetRealizable_;

    return realizable;
}

bool Foam::univariateMomentSet::isSubsetRealizable
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!realizabilityChecked_)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return subsetRealizable_;
}

bool Foam::univariateMomentSet::isOnMomentSpaceBoundary
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!realizabilityChecked_)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return onMomentSpaceBoundary_;
}

Foam::label Foam::univariateMomentSet::nRealizableMoments
(
    bool fatalErrorOnFailedRealizabilityTest
)
{
    if (!realizabilityChecked_)
    {
        checkRealizability(fatalErrorOnFailedRealizabilityTest);
    }

    return nRealizableMoments_;
}

Foam::scalar& Foam::univariateMomentSet::operator[](const label index)
{
    degenerate_ = false;
    fullyRealizable_ = true;
    subsetRealizable_ = true;
    onMomentSpaceBoundary_ = false;
    realizabilityChecked_ = false;

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
    return operator[](momentOrder);
}

Foam::scalar& Foam::univariateMomentSet::operator()(label momentOrder)
{
    return operator[](momentOrder);
}

// ************************************************************************* //
