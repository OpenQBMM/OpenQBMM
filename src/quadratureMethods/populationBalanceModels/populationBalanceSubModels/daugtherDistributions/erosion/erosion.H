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
    Foam::populationBalanceSubModels::daughterDistribution::erosion

Description
    Daughter distribution of fragments for erosion.

SourceFiles
    erosion.C

\*---------------------------------------------------------------------------*/

#ifndef erosion_H
#define erosion_H

#include "daughterDistribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace daughterDistributions
{

/*---------------------------------------------------------------------------*\
                    Class erosion Declaration
\*---------------------------------------------------------------------------*/

class erosion
:
    public daughterDistribution
{
private:
    
    // Private data

        //- Primary particle size
        dimensionedScalar primarySize_;


public:

        //- Runtime type information
        TypeName("erosion");


    // Constructors

        //- Construct from components
        erosion(const dictionary& dict);


    //- Destructor
    virtual ~erosion();


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
