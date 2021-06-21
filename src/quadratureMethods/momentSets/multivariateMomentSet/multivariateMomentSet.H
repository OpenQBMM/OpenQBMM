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
    Copyright (C) 2019 Alberto Passalacqua
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

Class
    Foam::multivariateMomentSet

Description
    Stores a multivariate moment set with the corresponding list of orders.

SourceFiles
    multivariateMomentSet.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateMomentSet_H
#define multivariateMomentSet_H

#include "momentSet.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multivariateMomentSet Declaration
\*---------------------------------------------------------------------------*/

class multivariateMomentSet
:
    public momentSet
{
public:

    // Constructors

        //- Construct from number of moments and initial value
        multivariateMomentSet
        (
            const label nMoments,
            const labelListList& momentOrders,
            const word& support,
            const scalar initValue = 0
        );

        //- Construct from scalarList
        multivariateMomentSet
        (
            const scalarList& m,
            const labelListList& momentOrders,
            const word& support
        );

    //- Destructor
    virtual ~multivariateMomentSet();


    // Member Functions

        //- Returns const access to the moment map
        inline const Map<label>& momentMap() const;

        //- Recalculate the moments from the quadrature nodes
        //void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multivariateMomentSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
