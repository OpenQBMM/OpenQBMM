/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2014 Bo Kong
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

    scalar     maxAbs_;

    //- Private member functions

    //- Disallow default bitwise copy construct
    hermiteQuadrature(const hermiteQuadrature&);

    //- Disallow default bitwise assignment
    void operator=(const hermiteQuadrature&);


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
    ) ;


    //- Return eigenvalues
    const scalarList& hermiteWeights() const
    {
        return origWei_;
    }

    //- Return eigenvectors
    const vectorList& hermiteAbscissas() const
    {
        return resAbs_;
    }

    vectorList& hermiteAbscissas()
    {
        return resAbs_;
    }

    const vectorList& hermiteOriginalAbscissas() const
    {
        return origAbs_;
    }

    const scalar& maxAbs() const
    {
        return maxAbs_;
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
