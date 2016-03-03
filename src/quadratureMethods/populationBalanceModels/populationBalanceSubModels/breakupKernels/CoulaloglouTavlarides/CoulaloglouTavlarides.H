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
    Foam::populationBalanceSubModels::breakupKernels::CoulaloglouTavlarides

Description
    Breakup kernel for liquid-liquid systems. 
    
    Reference
    \verbatim
        "Description of interaction processes in agitated liquid-liquid dis[ersions"
        C Coulaloglou, L Tavlarides
        Chemical Engineering Science 
        Volume 32, Pages 1289-1297, 1977
    \endverbatim 
    
SourceFiles
    CoulaloglouTavlarides.C

\*---------------------------------------------------------------------------*/

#ifndef CoulaloglouTavlarides_H
#define CoulaloglouTavlarides_H

#include "breakupKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{

/*---------------------------------------------------------------------------*\
                    Class CoulaloglouTavlarides Declaration
\*---------------------------------------------------------------------------*/

class CoulaloglouTavlarides
:
    public breakupKernel
{
private:
    
    // Private data

        //- Model coefficient
        const dimensionedScalar C1b_;

        //- Model coefficient
        const dimensionedScalar C2b_;
        
        //- Surface tension
        const dimensionedScalar surfaceT_;
        
    
public:

        //- Runtime type information
        TypeName("CoulaloglouTavlarides");


    // Constructors

        //- Construct from components
        CoulaloglouTavlarides(const dictionary& dict);


    //- Destructor
    virtual ~CoulaloglouTavlarides();


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
