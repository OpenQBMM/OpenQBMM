/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::daughterDistributions::symmetricFragmentation

Description
    Daughter distribution of fragments with symmetric fragmentation.

    The moment source term corresponding to the distribution is:

    \f[
        \overline{b}_k = 2^{(3-k)/3} \xi^k
    \f]

    \vartable
        k           |  Order of the moment [-]
        \xi         |  Particle size [m]
    \endvartable

SourceFiles
    symmetricFragmentation.C

\*---------------------------------------------------------------------------*/

#ifndef symmetricFragmentation_H
#define symmetricFragmentation_H

#include "daughterDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace daughterDistributions
{

/*---------------------------------------------------------------------------*\
                    Class symmetricFragmentation Declaration
\*---------------------------------------------------------------------------*/

class symmetricFragmentation
:
    public daughterDistribution
{
public:

        //- Runtime type information
        TypeName("symmetricFragmentation");


    // Constructors

        //- Construct from components
        symmetricFragmentation(const dictionary& dict);


    //- Destructor
    virtual ~symmetricFragmentation();


    // Member Functions

        //- Moment of the daughter distribution function
        virtual scalar mD
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
