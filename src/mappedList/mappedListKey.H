/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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
    along with OpenQBMM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mappedListKey

Description
    The index a mapped list finds its elements by: the orders, or indexes, of
    every element, and a table from them to the position of the element.

    An element is identified by its orders themselves. They are hashed to
    find where to look, and the orders stored there are compared with the
    ones asked for before the element is returned, so two different lists of
    orders are never taken for one another, whatever their values and however
    many dimensions they have.

    The keys used to pack the orders into one label, a decimal digit to a
    dimension, which made an order of ten or more land on the key of another
    list of orders, and limited a list to as many dimensions as a label has
    digits.

    Orders followed by zeros are the same element as the orders without
    them: (1) and (1 0 0) are the same element of a list of three
    dimensions, so an element can be asked for by its leading orders only,
    as moments(1) for the first moment along the first dimension.

    Where the orders of every dimension fit into 64 bits, as they do for any
    list of a few dimensions, the key is the orders packed, which is exact
    and is compared as one number; only a list of many dimensions falls back
    on the hash and the comparison of the orders.

    The orders and the table are held in flat lists, so an index, and a
    mapped list with it, is built and copied with a few allocations whatever
    its size: the inversions build and copy mapped lists in every cell, and
    the collision kernels look moments up for every pair of nodes.

SourceFiles
    mappedListKey.H

\*---------------------------------------------------------------------------*/

#ifndef mappedListKey_H
#define mappedListKey_H

#include "labelList.H"
#include "error.H"
#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mappedListKey Declaration
\*---------------------------------------------------------------------------*/

class mappedListKey
{
    // Private data

        //- Number of dimensions of the orders
        label nDimensions_;

        //- Number of elements
        label nElements_;

        //- Bits an order takes in a packed key, 0 when the orders of the
        //  list do not fit into one
        label bits_;

        //- Shift that takes a mixed key to a slot of the table
        label shift_;

        //- The orders of every element, nDimensions_ to an element, padded
        //  with zeros
        labelList orders_;

        //- Open-addressed table of positions, -1 where empty. Its size is a
        //  power of two at least twice the number of elements.
        labelList slots_;

        //- Key of the element in each slot of the table
        List<uint64_t> slotKeys_;


    // Private Member Functions

        //- Slot a key is first looked for at (Fibonacci hashing)
        label firstSlot(const uint64_t key) const
        {
            return label((key*11400714819323198485ULL) >> shift_);
        }

        //- Key of size orders from first, and whether they can be the orders
        //  of an element at all.
        //
        //  Where the orders fit, a dimension to bits_ bits, the key is the
        //  orders packed, which is exact: two different orders are two
        //  different keys, and an order too large for its bits belongs to
        //  no element. Otherwise the key is a hash of the orders (FNV-1a),
        //  which the orders themselves are compared with. Either way the
        //  orders are read to the number of dimensions of the list, so that
        //  orders followed by zeros are the orders without them.
        template<class Iterator>
        bool key(Iterator first, const label size, uint64_t& k) const
        {
            Iterator iter = first;

            if (bits_ > 0)
            {
                const uint64_t limit = uint64_t(1) << bits_;

                k = 0;

                for (label i = 0; i < size; ++i, ++iter)
                {
                    const uint64_t v = uint64_t(*iter);

                    if (i < nDimensions_)
                    {
                        if (v >= limit)
                        {
                            return false;
                        }

                        k |= v << (bits_*i);
                    }
                    else if (v != 0)
                    {
                        return false;
                    }
                }

                return true;
            }

            uint64_t h = 14695981039346656037ULL;

            const label n = min(size, nDimensions_);

            for (label i = 0; i < n; ++i, ++iter)
            {
                h = (h ^ uint64_t(*iter))*1099511628211ULL;
            }

            for (label i = n; i < nDimensions_; ++i)
            {
                h *= 1099511628211ULL;
            }

            k = h;

            return true;
        }

        //- Whether size orders from first are those of an element, the
        //  shorter of the two read with zeros past its end
        template<class Iterator>
        bool equal
        (
            const label elemi,
            Iterator first,
            const label size
        ) const
        {
            const label* stored = orders_.cdata() + elemi*nDimensions_;
            const label n = min(size, nDimensions_);

            Iterator iter = first;

            for (label i = 0; i < n; ++i, ++iter)
            {
                if (*iter != stored[i])
                {
                    return false;
                }
            }

            for (label i = n; i < nDimensions_; ++i)
            {
                if (stored[i] != 0)
                {
                    return false;
                }
            }

            for (label i = n; i < size; ++i, ++iter)
            {
                if (*iter != 0)
                {
                    return false;
                }
            }

            return true;
        }

        //- Position of the element of size orders from first, -1 if there
        //  is none
        template<class Iterator>
        label find(Iterator first, const label size) const
        {
            uint64_t k;

            if (!key(first, size, k))
            {
                return -1;
            }

            const label mask = slots_.size() - 1;

            for (label s = firstSlot(k); ; s = (s + 1) & mask)
            {
                const label elemi = slots_[s];

                if (elemi < 0)
                {
                    return -1;
                }

                if
                (
                    slotKeys_[s] == k
                 && (bits_ > 0 || equal(elemi, first, size))
                )
                {
                    return elemi;
                }
            }
        }

        //- Build the table from the orders of every element
        void build(const UList<labelList>& indexes)
        {
            nElements_ = indexes.size();
            nDimensions_ = 0;

            label maxOrder = 0;

            forAll(indexes, elemi)
            {
                nDimensions_ = max(nDimensions_, indexes[elemi].size());

                forAll(indexes[elemi], dimi)
                {
                    if (indexes[elemi][dimi] < 0)
                    {
                        FatalErrorInFunction
                            << "Element " << elemi << " has the negative "
                            << "orders " << indexes[elemi] << "." << nl
                            << exit(FatalError);
                    }

                    maxOrder = max(maxOrder, indexes[elemi][dimi]);
                }
            }

            // The bits the largest order takes, and whether the orders of
            // every dimension fit into one key of 64 bits
            bits_ = 1;

            while ((label(1) << bits_) <= maxOrder)
            {
                bits_++;
            }

            if (bits_*nDimensions_ > 64)
            {
                bits_ = 0;
            }

            orders_.resize_nocopy(nElements_*nDimensions_);
            orders_ = Zero;

            forAll(indexes, elemi)
            {
                forAll(indexes[elemi], dimi)
                {
                    orders_[elemi*nDimensions_ + dimi] = indexes[elemi][dimi];
                }
            }

            label nSlots = 2;
            label logSlots = 1;

            while (nSlots < 2*nElements_)
            {
                nSlots *= 2;
                logSlots++;
            }

            shift_ = 64 - logSlots;

            slots_.resize_nocopy(nSlots);
            slots_ = -1;
            slotKeys_.resize_nocopy(nSlots);
            slotKeys_ = 0;

            const label mask = nSlots - 1;

            forAll(indexes, elemi)
            {
                const labelList& orders = indexes[elemi];

                uint64_t k;
                key(orders.cbegin(), orders.size(), k);

                label s = firstSlot(k);

                while (slots_[s] >= 0)
                {
                    const label other = slots_[s];

                    if
                    (
                        slotKeys_[s] == k
                     && (
                            bits_ > 0
                         || equal(other, orders.cbegin(), orders.size())
                        )
                    )
                    {
                        FatalErrorInFunction
                            << "Element " << elemi << " of orders " << orders
                            << " has the orders of element " << other << "."
                            << nl << exit(FatalError);
                    }

                    s = (s + 1) & mask;
                }

                slots_[s] = elemi;
                slotKeys_[s] = k;
            }
        }


public:

    // Constructors

        //- Construct empty
        mappedListKey()
        {
            build(UList<labelList>());
        }

        //- Construct from the orders of every element, stopping if two
        //  elements have the same orders
        explicit mappedListKey(const UList<labelList>& indexes)
        {
            build(indexes);
        }


    // Member Functions

        //- Number of dimensions of the orders
        label nDimensions() const noexcept
        {
            return nDimensions_;
        }

        //- Number of elements
        label size() const noexcept
        {
            return nElements_;
        }

        //- Position of the element of a list of orders, -1 if there is none
        label find(const UList<label>& orders) const
        {
            return find(orders.cbegin(), orders.size());
        }

        //- Position of the element of orders given one by one, -1 if there
        //  is none
        label find(std::initializer_list<label> orders) const
        {
            return find(orders.begin(), label(orders.size()));
        }

        //- The orders of every element, padded with zeros to the number of
        //  dimensions
        labelListList indexes() const
        {
            labelListList result(nElements_);

            forAll(result, elemi)
            {
                result[elemi] =
                    SubList<label>
                    (
                        orders_,
                        nDimensions_,
                        elemi*nDimensions_
                    );
            }

            return result;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
