/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019 Alberto Passalacqua
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
    Foam::PDFTransportModel

Description
    Abstract class for the run-time selection of the PDF transport model.

SourceFiles
    PDFTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef PDFTransportModel_H
#define PDFTransportModel_H

#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PDFTransportModel Declaration
\*---------------------------------------------------------------------------*/

class PDFTransportModel
{
    // Private data

        //- Name of the PDFTransportModel
        const word name_;


protected:

    // Protected data

        //- Constant reference to the mesh
        const fvMesh& mesh_;


public:

    // Constructors

        //- Construct from components
        PDFTransportModel
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construct
        PDFTransportModel(const PDFTransportModel&) = delete;


    //- Destructor
    virtual ~PDFTransportModel();


    // Member Functions

        //- Solve PDF transport equation
        virtual void solve() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PDFTransportModel&) = delete;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
