/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 Alberto Passalacqua
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

#include "lognormalEQMOM.H"
#include "scalar.H"
#include "scalarMatrices.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lognormalEQMOM, 0);

    addToRunTimeSelectionTable
    (
        extendedMomentInversion, 
        lognormalEQMOM, 
        dictionary
    ); 
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lognormalEQMOM::lognormalEQMOM
(
    const dictionary& dict,
    const univariateMomentSet& moments,
    const label nSecondaryNodes
)
:
    extendedMomentInversion(dict, moments, nSecondaryNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lognormalEQMOM::~lognormalEQMOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::lognormalEQMOM::secondaryAbscissa
(
    scalar primaryAbscissa,
    scalar secondaryAbscissa,
    scalar sigma
)
{
    return primaryAbscissa*exp(Foam::sqrt(2.0)*secondaryAbscissa*sigma);
}

void Foam::lognormalEQMOM::momentsStarToMoments
(
    scalar sigma,
    univariateMomentSet& moments
)
{
    scalar z = exp(sqr(sigma)/2.0);

    forAll(moments, momentI)
    {
        moments[momentI] = momentsStar_[momentI]*pow(z, momentI*momentI);
    }
}

void Foam::lognormalEQMOM::momentsToMomentsStar(scalar sigma)
{
    scalar z = exp(-sqr(sigma)/2.0);

    forAll(moments_, momentI)
    {
        momentsStar_[momentI] = moments_[momentI]*pow(z, momentI*momentI);
    }
}

void Foam::lognormalEQMOM::recurrenceRelation
(
    scalarDiagonalMatrix& a, 
    scalarDiagonalMatrix& b,
    scalar primaryAbscissa
)
{
//  Unnecessary re-initialization because we set a to be full of zeros and b is
//  recalculated
//     forAll(a, aI)
//     {
//         a[aI] = 0.0;
//     }
//     
//     forAll{b, bI)
//     {
//         b[bI] = 0.0;  
//     }

    forAll(b, bI)
    {
        if (bI == 0)
        {
            b[bI] = sqrt(Foam::constant::mathematical::pi); //tgamma(0.5);
        }
        else 
        {
            b[bI] = scalar(bI)/2.0;
        }
    }
}

Foam::scalar Foam::lognormalEQMOM::sigmaMax()
{
    scalar sigmaZeta1 
        = sqrt(2.0*log(sqrt(moments_[0]*moments_[2]/(sqr(moments_[0])))));
    
    scalar sigmaZeta2 
            = sqrt(2.0*log(sqrt(moments_[1]*moments_[3]/(sqr(moments_[2])))));

    return min(sigmaZeta1, sigmaZeta2);
}

// ************************************************************************* //
