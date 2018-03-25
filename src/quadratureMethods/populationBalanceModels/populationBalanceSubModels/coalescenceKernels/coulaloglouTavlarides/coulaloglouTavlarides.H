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
    Foam::populationBalanceSubModels::coalescenceKernels::coulaloglouTavlarides

Description
    coulaloglouTavlarides kernel supplemented kernel for 
    
    References
    \verbatim
        "Simulation of polydisperse multiphase systems using population balance
        and example application to bubbly flows"
        Antonio Buffo, Daniele Marchisio, Marco Vanni, Peter Renze
        Chemical Engineering Research and Design
        Volume 91, Pages 1859-1875, 2013
    \endverbatim 
    
    \verbatim
        "Description of interaction processes in agitated liquid-liquid
        dispersions"
        C A Coulaloglou, L L Tavlrides
        Chemical Engineering Science
        Volume 32, Pages 1289-1297, 1977
    \endverbatim 

SourceFiles
    coulaloglouTavlarides.C

\*---------------------------------------------------------------------------*/

#ifndef coulaloglouTavlarides_H
#define coulaloglouTavlarides_H

#include "coalescenceKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace coalescenceKernels
{

/*---------------------------------------------------------------------------*\
                    Class coulaloglouTavlarides Declaration
\*---------------------------------------------------------------------------*/

class coulaloglouTavlarides
:
    public coalescenceKernel
{
    // Private data
        //- Coefficients
        scalar etaC_;
public:

        //- Runtime type information
        TypeName("coulaloglouTavlarides");


    // Constructors

        //- Construct from components
        coulaloglouTavlarides(const dictionary& dict);


    //- Destructor
    virtual ~coulaloglouTavlarides();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<volScalarField> Ka
        (
            const volScalarField& abscissa1,
            const volScalarField& abscissa2
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coalescenceKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
