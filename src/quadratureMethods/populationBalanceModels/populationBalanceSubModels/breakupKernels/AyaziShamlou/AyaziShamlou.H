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
    Foam::populationBalanceSubModels::breakupKernels::AyaziShamlou

Description
    Breakup kernel for liquid-solid systems.
    
    Reference
    \verbatim
        "Growth-independent breakage frequency of protein precipitates in 
        turbulently agitated bioreactors"
        P Ayazi Shamlou, S Stravrinides, N Titchener-Hooker, M Hoare
        Chemical Engineering Science
        Volume 49, Pages 2647-2656, 1994
    \endverbatim 

SourceFiles
    AyaziShamlou.C

\*---------------------------------------------------------------------------*/

#ifndef AyaziShamlou_H
#define AyaziShamlou_H

#include "breakupKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{

/*---------------------------------------------------------------------------*\
                    Class AyaziShamlou Declaration
\*---------------------------------------------------------------------------*/

class AyaziShamlou
:
    public breakupKernel
{
private:
    
    // Private data

        //- Hamaker constant
        dimensionedScalar A_;
        
        //- Fractal dimension of the aggregates
        dimensionedScalar df_;
        
        //- Distance between primary particles
        dimensionedScalar H0_;        
        
        //- Primary particle size
        dimensionedScalar primarySize_;
    
public:

        //- Runtime type information
        TypeName("AyaziShamlou");


    // Constructors

        //- Construct from components
        AyaziShamlou(const dictionary& dict);


    //- Destructor
    virtual ~AyaziShamlou();


    // Member Functions

        //- Breakup kernel
        virtual tmp<volScalarField> Kb
        (
            const volScalarField& abscissa
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
