/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
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
    Foam::populationBalanceSubModels::daughterDistributions::oneQuarterMassRatio

Description
    Daughter distribution of fragments with one quarter mass ratio
    fragmentation.

    The moment source term corresponding to the distribution is:

    \f[
        \overline{b}_k = \frac{4^{k/3} + 1}{5^{k/3}} \xi^k
    \f]

    where

    \vartable
        k           |  Order of the moment [-]
        \xi         |  Particle size [m]
    \endvartable

SourceFiles
    oneQuarterMassRatio.C

\*---------------------------------------------------------------------------*/

#ifndef oneQuarterMassRatio_H
#define oneQuarterMassRatio_H

#include "daughterDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace daughterDistributions
{

/*---------------------------------------------------------------------------*\
                    Class oneQuarterMassRatio Declaration
\*---------------------------------------------------------------------------*/

class oneQuarterMassRatio
:
    public daughterDistribution
{
public:

        //- Runtime type information
        TypeName("oneQuarterMassRatio");


    // Constructors

        //- Construct from components
        oneQuarterMassRatio(const dictionary& dict);


    //- Destructor
    virtual ~oneQuarterMassRatio();


    // Member Functions

        //- Moment of the daughter distribution function (length based)
        virtual scalar mD
        (
            const label& order,
            const scalar& abscissa
        ) const;

        //- Moment of the daughter distribution function (mass/volume based)
        virtual scalar mDMass
        (
            const label& order,
            const scalar& abscissa
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterDistributions
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
