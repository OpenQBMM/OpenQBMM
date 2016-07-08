/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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
    Foam::populationBalanceModels::daughterDistribution::fullFragmentation

Description
    Daughter distribution of fragments with full fragmentation.

SourceFiles
    fullFragmentation.C

\*---------------------------------------------------------------------------*/

#ifndef fullFragmentation_H
#define fullFragmentation_H

#include "daughterDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace daughterDistributions
{

/*---------------------------------------------------------------------------*\
                    Class fullFragmentation Declaration
\*---------------------------------------------------------------------------*/

class fullFragmentation
:
    public daughterDistribution
{
private:

    // Private data

        //- Primary particle size
        dimensionedScalar primarySize_;


public:

        //- Runtime type information
        TypeName("fullFragmentation");


    // Constructors

        //- Construct from components
        fullFragmentation(const dictionary& dict);


    //- Destructor
    virtual ~fullFragmentation();


    // Member Functions

        //- Moment of the daughter distribution function
        virtual tmp<volScalarField> mD
        (
            const label order,
            const volScalarField& abscissa
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace daughterDistributions
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
