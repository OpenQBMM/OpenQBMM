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

#include "gammaEQMOM.H"
#include "scalar.H"
#include "scalarMatrices.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gammaEQMOM, 0);

    addToRunTimeSelectionTable
    (
        extendedMomentInversion, 
        gammaEQMOM, 
        dictionary
    ); 
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gammaEQMOM::gammaEQMOM
(
    const dictionary& dict,
    const label nMoments,    
    const label nSecondaryNodes
)
:
    extendedMomentInversion(dict, nMoments, nSecondaryNodes)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gammaEQMOM::~gammaEQMOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::gammaEQMOM::secondaryAbscissa
(
    scalar primaryAbscissa,
    scalar secondaryAbscissa,
    scalar sigma
)
{
    return sigma*secondaryAbscissa;
}

void Foam::gammaEQMOM::momentsStarToMoments
(
    scalar sigma,
    univariateMomentSet& moments,
    const univariateMomentSet& momentsStar
)
{
    label nMom = moments.size();
    
    if (nMom >= 8)
    {
        FatalErrorIn
        (
            "Foam::gammaEQMOM::momentsStarToMoments\n"
            "(\n"
            "   scalar sigma,\n"
            "   univariateMomentSet& moments,\n"
            "   const univariateMomentSet& momentsStar\n"
            ")"
        )   << "Moment transformation not implemented."
            << abort(FatalError);
    }
    
    moments[0] = momentsStar[0];
    moments[1] = momentsStar[1];
    moments[2] = momentsStar[2] + sigma*momentsStar[1];
    
    if (nMom >= 5)
    {
        moments[3] = momentsStar[3] + 3.0*sigma*momentsStar[2] 
                + 2.0*sqr(sigma)*momentsStar[1];

        moments[4] = momentsStar[4] + 6.0*sigma*momentsStar[3]
                + 11.0*sqr(sigma)*momentsStar[2] 
                + 6.0*pow3(sigma)*momentsStar[1];
    }
    
    if (nMom >= 7)
    {       
        moments[5] = momentsStar[5] + 10.0*sigma*momentsStar[4]
                + 35.0*sqr(sigma)*momentsStar[3]
                + 50.0*pow3(sigma)*momentsStar[2]
                + 24.0*pow4(sigma)*momentsStar[1];

        moments[6] = momentsStar[6] + 15.0*sigma*momentsStar[5]
                + 85.0*sqr(sigma)*momentsStar[4]
                + 225.0*pow3(sigma)*momentsStar[3]
                + 274.0*pow4(sigma)*momentsStar[2]
                + 120.0*pow5(sigma)*momentsStar[1];
    }
    
    if (nMom >= 8)
    {
        FatalErrorIn
        (
            "Foam::gammaEQMOM::momentsStarToMoments\n"
            "(\n"
            "   scalar sigma,\n"
            "   univariateMomentSet& moments,\n"
            "   const univariateMomentSet& momentsStar\n"
            ")"
        )   << "Moment transformation not implemented."
            << abort(FatalError);
    }
}

void Foam::gammaEQMOM::momentsToMomentsStar
(
    scalar sigma,
    const univariateMomentSet& moments,
    univariateMomentSet& momentsStar
)
{
    label nMom = moments.size();

    if (nMom >= 8)
    {
        FatalErrorIn
        (
            "Foam::gammaEQMOM::momentsToMomentsStar\n"
            "(\n"
            "   scalar sigma,\n"
            "   const univariateMomentSet& moments,\n"
            "   univariateMomentSet& momentsStar\n"
            ")"
        )   << "Moment transformation not implemented."
            << abort(FatalError);
    }
    
    momentsStar[0] = moments[0];
    momentsStar[1] = moments[1];
    momentsStar[2] = momentsStar[2] - sigma*momentsStar[1];    

    if (nMom >= 5)
    {
        momentsStar[3] = moments[3] - 2.0*sigma*moments[2];

        momentsStar[4] = moments[4] - 6.0*sigma*moments[3]
                + sqr(sigma)*moments[2] + 5.0*pow3(sigma)*moments[1];
    }
    
    if (nMom >= 7)
    {
        momentsStar[5] = moments[5] - 10.0*sigma*moments[4]
                + 25.0*sqr(sigma)*moments[3] + 10.0*pow3(sigma)*moments[2]
                - 24.0*pow4(sigma)*moments[1];
                
       momentsStar[6] = moments[6] - 15.0*sigma*moments[5]
                + 65.0*sqr(sigma)*moments[4] - 90.0*pow3(sigma)*moments[3]
                - 59.0*pow4(sigma)*moments[2] + 89.0*pow5(sigma)*moments[1];
    }
}

Foam::scalar Foam::gammaEQMOM::m2N
(
    scalar sigma, 
    univariateMomentSet momentsStar
)
{   
    label nMomentsStar = momentsStar.size();
    
    if (momentsStar.nRealizableMoments() >= nMomentsStar - 1)
    {
        univariateMomentSet m(nMomentsStar, 0.0);
        momentsStarToMoments(sigma, m, momentsStar);
        
        return m.last();
    }
    
    return GREAT;
}

void Foam::gammaEQMOM::recurrenceRelation
(
    scalarDiagonalMatrix& a, 
    scalarDiagonalMatrix& b,
    scalar primaryAbscissa,
    scalar sigma
)
{
    scalar alpha = primaryAbscissa/sigma - 1.0;
    
    forAll(a, aI)
    {
        a[aI] = (2.0*scalar(aI) + alpha + 1.0);
    }

    b[0] = gamma(1.0 + alpha);

    for (label bI = 1; bI < b.size(); bI++)
    {
        b[bI] = scalar(bI)*(scalar(bI) + alpha);
    }
}

Foam::scalar Foam::gammaEQMOM::sigmaMax(const univariateMomentSet& moments)
{  
    scalar sigmaZeta1 = 
        (moments[0]*moments[2] - moments[1]*moments[1])/(moments[0]*moments[1]);
    
    return sigmaZeta1;//min(sigmaZeta1, sigmaZeta2);
}

// ************************************************************************* //
