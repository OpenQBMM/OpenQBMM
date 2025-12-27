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

#include "momentSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentSet::momentSet
(
    const label nMoments,
    const label nDimensions,
    const labelListList& momentOrders,
    const List<supportType>& supports,
    const scalar smallM0,
    const scalar smallZeta,
    const scalar initialValue
)
:
    mappedList(nMoments, momentOrders, initialValue),
    nDimensions_(nDimensions),
    momentOrders_(momentOrders),
    supports_(supports),
    smallM0_(smallM0),
    smallZeta_(smallZeta)
{
    validate();
}

Foam::momentSet::momentSet
(
    const scalarList& m,
    const label nDimensions,
    const labelListList& momentOrders,
    const List<supportType>& supports,
    const scalar smallM0,
    const scalar smallZeta
)
:
    mappedList(m, momentOrders),
    nDimensions_(nDimensions),
    momentOrders_(momentOrders),
    supports_(supports),
    smallM0_(smallM0),
    smallZeta_(smallZeta)
{
    validate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentSet::~momentSet()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::momentSet::setSize
(
    const label newSize,
    labelListList& newMomentOrders
)
{
    // Check if new size is valid
    if (newSize < 2)
    {
        FatalErrorInFunction
            << "The new size of the moment set must be greater than 1." << endl
            << abort(FatalError);
    }

    // Do not resize if the size is unchanged
    if (newSize == nMoments())
    {
        return;
    }

    Foam::mappedList<scalar>::setSize(newSize, newMomentOrders);

    // Check that the new moment orders are valid
    if (newMomentOrders.size() != newSize)
    {
        FatalErrorInFunction
            << "The size of the new moment orders list is inconsistent "<< endl
            << "with the new size of the moment set." << endl
            << "    New size of moment set: " << newSize << endl
            << "    Size of new moment orders list: " << newMomentOrders.size()
            << endl
            << abort(FatalError);
    }

    momentOrders_ = newMomentOrders;
    validateMomentOrders();
}

void Foam::momentSet::resize
(
    const label newSize,
    labelListList& newMomentOrders
)
{
    (*this).setSize(newSize, newMomentOrders);
}

void Foam::momentSet::validateSupports() const
{
    forAll(supports_, i)
    {
        if (supports_[i] != supportType::R
         && supports_[i] != supportType::RPlus
         && supports_[i] != supportType::ZeroOne)
        {
            FatalErrorInFunction
                << "The specified support for dimension " << i
                << " is invalid." << nl
                << "    Valid supports are: R, RPlus and 01." << nl
                << "    Moment set: " << (*this)
                << abort(FatalError);
        }
    }
}

void Foam::momentSet::validateMomentOrders() const
{
    forAll(momentOrders_, i)
    {
        if (momentOrders_[i].size() != nDimensions_)
        {
            FatalErrorInFunction
                << "The moment order " << momentOrders_[i]
                << " does not have the correct number of dimensions." << nl
                << "    Expected number of dimensions: " << nDimensions_ << nl
                << "    Moment set: " << (*this)
                << abort(FatalError);
        }
    }
}

void Foam::momentSet::validate() const
{
    validateSupports();
    validateMomentOrders();

    // Ensure that moments are finite and throw an exception otherwise
    forAll((*this), i)
    {
        if (!std::isfinite((*this)[i]))
        {
            FatalErrorInFunction
                << "Moment " << i << " is not finite: " << (*this)[i] << nl
                << "    Moment set: " << (*this)
                << abort(FatalError);
        }
    }
}

// ************************************************************************* //
