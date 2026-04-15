/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "fixedFaceFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedFaceFvPatchScalarField::fixedFaceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF),
    fixProc_(false)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;

    if (UPstream::parRun())
    {
        #if (OPENFOAM >= 2512)
        fixProc_ =
        (
            UPstream::find_first(p.size() > 0, UPstream::worldComm)
         == UPstream::myProcNo(UPstream::worldComm)
        );
        #else
        label minProc = p.size() > 0 ? Pstream::myProcNo() : Pstream::nProcs();

        Foam::reduce(minProc, minOp<label>());

        fixProc_ = (minProc == Pstream::myProcNo());
        #endif
    }
    else
    {
        fixProc_ = (p.size() > 0);
    }

    if (fixProc_)
    {
        valueFraction()[0] = 1.0;
    }
}


Foam::fixedFaceFvPatchScalarField::fixedFaceFvPatchScalarField
(
    const this_bctype& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    fixProc_(ptf.fixProc_)
{}


Foam::fixedFaceFvPatchScalarField::fixedFaceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
    fixProc_(false)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;

    if (UPstream::parRun())
    {
        #if (OPENFOAM >= 2512)
        fixProc_ =
        (
            UPstream::find_first(p.size() > 0, UPstream::worldComm)
         == UPstream::myProcNo(UPstream::worldComm)
        );
        #else
        label minProc = p.size() > 0 ? Pstream::myProcNo() : Pstream::nProcs();

        Foam::reduce(minProc, minOp<label>());

        fixProc_ = (minProc == Pstream::myProcNo());
        #endif
    }
    else
    {
        fixProc_ = (p.size() > 0);
    }

    if (fixProc_)
    {
        valueFraction()[0] = 1.0;
    }
}


Foam::fixedFaceFvPatchScalarField::fixedFaceFvPatchScalarField
(
    const this_bctype& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    fixProc_(ptf.fixProc_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFaceFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (fixProc_)
    {
        valueFraction()[0] = 1.0;
    }

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::fixedFaceFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
    //writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFaceFvPatchScalarField
    );
}

// ************************************************************************* //
