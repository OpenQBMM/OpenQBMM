/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/

#include "univariateMomentInversion.H"
#include "IOmanip.H"

#include "eigenSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(univariateMomentInversion, 0);
    defineRunTimeSelectionTable(univariateMomentInversion, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::univariateMomentInversion::univariateMomentInversion
(
    const dictionary& dict
)
:
    nNodes_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::univariateMomentInversion::~univariateMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::univariateMomentInversion::JacobiMatrix
(
    univariateMomentSet& moments,
    scalarSquareMatrix& z,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    correctRecurrence(moments, minKnownAbscissa, maxKnownAbscissa);

    for (label i = 0; i < nNodes_ - 1; i++)
    {
        z[i][i] = moments.alphaRecurrence()[i];
        z[i][i+1] = Foam::sqrt(moments.betaRecurrence()[i+1]);
        z[i+1][i] = z[i][i+1];
    }

    z[nNodes_ - 1][nNodes_ - 1] = moments.alphaRecurrence()[nNodes_ - 1];
}

void Foam::univariateMomentInversion::invert
(
    univariateMomentSet& moments,
    scalarList& weights,
    scalarList& abscissae,
    const scalar minKnownAbscissa,
    const scalar maxKnownAbscissa
)
{
    if (moments.isDegenerate())
    {
        return;
    }

    if (moments[0] < SMALL)
    {
        nNodes_ = 0;

        weights.setSize(nNodes_);
        abscissae.setSize(nNodes_);

        return;
    }

    calcNQuadratureNodes(moments, weights, abscissae);

    if (nInvertibleMoments_ == 2)
    {
        weights[0] = moments[0];
        abscissae[0] = moments[1]/moments[0];

        return;
    }

    scalarSquareMatrix z(nNodes_, scalar(0));
    JacobiMatrix(moments, z, minKnownAbscissa, maxKnownAbscissa);

    // Computing weights and abscissae
    eigenSolver zEig(z, true);

    // Computing weights and abscissae
    for (label i = 0; i < nNodes_; i++)
    {
        weights[i] = moments[0]*sqr(zEig.eigenvectors()[0][i]);
        abscissae[i] = zEig.eigenvaluesRe()[i];
    }
}


// ************************************************************************* //
