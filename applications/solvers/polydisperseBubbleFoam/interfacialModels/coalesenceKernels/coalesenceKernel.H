/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 Alberto Passalacqua
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
    Foam::populationBalanceSubModels::coalesenceKernel

Description
    Abstract class for coalesence kernels.

SourceFiles
    coalesenceKernel.C

\*---------------------------------------------------------------------------*/

#ifndef coalesenceKernel_H
#define coalesenceKernel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "coalesenceFrequencyKernel.H"
#include "coalesenceEfficiencyKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{

/*---------------------------------------------------------------------------*\
                    Class coalesenceKernel Declaration
\*---------------------------------------------------------------------------*/

class coalesenceKernel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        coalesenceKernel(const coalesenceKernel&);

        //- Disallow default bitwise assignment
        void operator=(const coalesenceKernel&);


    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Coefficient of coalesence kernel
        const dimensionedScalar Ca_;

        //- Coalesence frequency model
        autoPtr<coalesenceFrequencyKernel> frequency_;

        //- Coalesence effeciency model
        autoPtr<coalesenceEfficiencyKernel> efficiency_;


public:

    // Constructor
    coalesenceKernel
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~coalesenceKernel();


    // Member Functions

        //- Coalesence kernel
        tmp<volScalarField> Ka
        (
            const volScalarField& abscissa1,
            const volScalarField& abscissa2
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
