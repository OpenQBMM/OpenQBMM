/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Class
    Foam::phasePairKey

Description
    An ordered or unorder pair of phase names.
    Typically specified as follows.
    \verbatim
       (phase1 and phase2)  // unordered
       (phase1 in  phase2)  // ordered
    \endverbatim

SourceFiles
    phasePairKey.C

\*---------------------------------------------------------------------------*/

#ifndef phasePairKey_H
#define phasePairKey_H

#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class phasePairKey;

bool operator==(const phasePairKey& a, const phasePairKey& b);
bool operator!=(const phasePairKey& a, const phasePairKey& b);

Istream& operator>>(Istream& is, phasePairKey& key);
Ostream& operator<<(Ostream& os, const phasePairKey& key);

/*---------------------------------------------------------------------------*\
                         Class phasePairKey Declaration
\*---------------------------------------------------------------------------*/

class phasePairKey
:
    public Pair<word>
{
    // Private Data

        //- Flag to indicate whether ordering is important
        bool ordered_;


public:

    // Constructors

        //- Default construct
        phasePairKey() = default;

        //- Construct from names and (optional) ordering flag
        phasePairKey
        (
            const word& name1,
            const word& name2,
            const bool ordered = false
        );


    //- Destructor
    virtual ~phasePairKey() = default;


    // Member Functions

        //- Return the ordered flag
        bool ordered() const noexcept
        {
            return ordered_;
        }


    // Friend Operators

        //- Test for equality
        friend bool operator==(const phasePairKey& a, const phasePairKey& b);

        //- Test for inequality
        friend bool operator!=(const phasePairKey& a, const phasePairKey& b);

        //- Read from Istream
        friend Istream& operator>>(Istream& is, phasePairKey& key);

        //- Write to Ostream
        friend Ostream& operator<<(Ostream& os, const phasePairKey& key);


    // Hashing

        //- Hashing functor for phasePairKey
        struct hasher
        {
            unsigned operator()(const phasePairKey& key) const
            {
                string::hash op;
                if (key.ordered())
                {
                    // return Pair<word>::hasher()(key);
                    return op(key.first(), op(key.second()));
                }
                else
                {
                    // return Pair<word>::symmHasher()(key);
                    // -> almost the same
                    return unsigned(op(key.first()) + op(key.second()));
                }
            }
        };

        //- Alternative name for functor
        using hash = phasePairKey::hasher;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hashing for phasePairKey
template<> struct Hash<phasePairKey> : phasePairKey::hasher {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
