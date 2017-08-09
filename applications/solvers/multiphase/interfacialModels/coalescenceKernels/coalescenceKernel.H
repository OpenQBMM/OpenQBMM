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
    Foam::populationBalanceSubModels::coalescenceKernel

Description
    Class to calculate the coalescence rate for bubbles based on coalescence
    frequency and efficiency.

SourceFiles
    coalescenceKernel.C

\*---------------------------------------------------------------------------*/

#ifndef coalescenceKernel_H
#define coalescenceKernel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "coalescenceFrequencyKernel.H"
#include "coalescenceEfficiencyKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class coalescenceKernel Declaration
\*---------------------------------------------------------------------------*/

class coalescenceKernel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        coalescenceKernel(const coalescenceKernel&);

        //- Disallow default bitwise assignment
        void operator=(const coalescenceKernel&);


    // Private data

        //- Dictionary
        const dictionary& dict_;

        //- Mesh
        const fvMesh& mesh_;

        //- Coefficient of coalescence kernel
        const dimensionedScalar Ca_;

        //- Coalescence frequency model
        autoPtr<coalescenceFrequencyKernel> frequency_;

        //- Coalescence effeciency model
        autoPtr<coalescenceEfficiencyKernel> efficiency_;


public:

    // Constructor
    coalescenceKernel
    (
        const dictionary& dict,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~coalescenceKernel();


    // Member Functions

        //- Coalescence kernel
        tmp<volScalarField> Ka
        (
            const label nodei,
            const label nodej
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
