/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alberto Passalacqua
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
    Fokker-Planck mixing model.

Description
    Brownian aggregation kernel supplemented kernel for particles in turbulent
    flows in conditions when the particle size is smaller than the local
    Kolmogorov scale.

    References
    \verbatim
        "Computational Models for Turbulent Reacting Flows"
        R O Fox
        Cambridge University Press
    \endverbatim
    
    \verbatim
        "Solution of the Fokker-Planck Equation Using the Extended
        Quadrature Method of Moments"
        E Madadi-Kandjani and A Passalacqua 
        2015 AIChE Annual Meeting, Salt Lake City, UT
    \endverbatim

SourceFiles
    FokkerPlanck.C

\*---------------------------------------------------------------------------*/

#ifndef FokkerPlanck_H
#define FokkerPlanck_H

#include "mixingKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace mixingKernels
{

/*---------------------------------------------------------------------------*\
                           Class FokkerPlanck Declaration
\*---------------------------------------------------------------------------*/

class FokkerPlanck
:
    public mixingKernel
{
public:

        //- Runtime type information
        TypeName("FokkerPlanck");


    // Constructors

        //- Construct from components
        FokkerPlanck(const dictionary& dict);


    //- Destructor
    virtual ~FokkerPlanck();


    // Member Functions

        //- Aggregation kernel
        virtual tmp<fvScalarMatrix> K
        (
            const volUnivariateMoment& moment,
            const volUnivariateMomentFieldSet& moments
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixingKernels
} // End namespace mixingSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
