/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2010-2014 by Bo Kong and 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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
    Foam::hermiteQuadrature

Description
   Use spectral methods to generate samples for multivaraites Gaussian
   distribution

\*---------------------------------------------------------------------------*/

#ifndef hermiteQuadrature_H
#define hermiteQuadrature_H

#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "scalarList.H"
#include "vectorList.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class hermiteQuadrature Declaration
\*---------------------------------------------------------------------------*/

class hermiteQuadrature
{
    // Private data

        static const scalar thetaLimit_;

        //- Number of weights and abscissas
        const label nDim_;

        const label nOrder_;

        const label nTolNodes_;

        scalarList origWei_;

        vectorList origAbs_;

        vectorList resAbs_;

        scalar maxAbs_;


public:

    // Constructors

        //- Construct from a scalarDiagonalMatrix
        hermiteQuadrature
        (
            const label& nDim,
            const label& nOrder
        );


        void calcHermiteQuadrature
        (
            const vector& mu,
            const symmTensor& Pp
        );

        //- Disallow default bitwise copy construct
        hermiteQuadrature(const hermiteQuadrature&) = delete;

    // Access Methods

        //- Return eigenvalues
        const scalarList& hermiteWeights() const
        {
            return origWei_;
        }

        //- Return eigenvectors
        const vectorList& hermiteAbscissae() const
        {
            return resAbs_;
        }

        vectorList& hermiteAbscissae()
        {
            return resAbs_;
        }

        const vectorList& hermiteOriginalAbscissae() const
        {
            return origAbs_;
        }

        const scalar& maxAbs() const
        {
            return maxAbs_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hermiteQuadrature&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
