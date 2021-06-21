/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2020 Alberto Passalacqua - Implemented Cunningham correction.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::dragModels::CunninghamSchillerNaumann

Description
    Schiller and Naumann drag model for isolated particles with Cunnigham
    correction factor.

    References
    \verbatim
        "A drag coefficient correlation"
        Schiller, L., Naumann Z.,
        Zeitschrift des Vereins Deutscher Ingenieure
        Volume 77, Page 318-320, 1935
    \endverbatim

    \verbatim
        "On the velocity of steady fall of spherical particles through fluid 
        medium"
        Cunningham, E., Larmor, J.
        Proceedings of the Royal Society of London. Series A
        Volume 83, Pages 357–365, 1910
    \endverbatim

    \verbatim
        "Definitive equations for the fluid resistance of spheres"
        Davies, C.N.
        Proceedings of the Physical Society
        Volume 57, Pages 259–270, 1945
    \endverbatim

SourceFiles
    CunninghamSchillerNaumann.C

\*---------------------------------------------------------------------------*/

#ifndef CunninghamSchillerNaumann_H
#define CunninghamSchillerNaumann_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

class SchillerNaumann;

/*---------------------------------------------------------------------------*\
                Class CunninghamSchillerNaumann Declaration
\*---------------------------------------------------------------------------*/

class CunninghamSchillerNaumann : public dragModel
{
    // Private data

        //- Ergun drag model
        autoPtr<SchillerNaumann> SchillerNaumann_;

        //- Residual Reynolds Number
        const dimensionedScalar residualRe_;

        //- Cunningham coefficient A1
        const dimensionedScalar A1_;

        //- Cunningham coefficient A2
        const dimensionedScalar A2_;

        //- Cunningham coefficient A3
        const dimensionedScalar A3_;

        //- Gas molecular weight. Defaults to air.
        const dimensionedScalar M_;


public:

    //- Runtime type information
    TypeName("CunninghamSchillerNaumann");


    // Constructors

        //- Construct from a dictionary and a phase pair
        CunninghamSchillerNaumann
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~CunninghamSchillerNaumann();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe(const label, const label) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
