/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the total number of moments
Foam::label Foam::conditionalMomentInversion::nMoments()
{
    return nMoments_;
}

//- Return the list of nNodes
const Foam::labelList& Foam::conditionalMomentInversion::nNodes()
{
    return nNodes_;
}

//- Return mapped list of weights
const Foam::PtrList<Foam::mappedList< Foam::scalar>>&
Foam::conditionalMomentInversion::weights() const
{
    return weights_;
}

//- Return mapped list of abscissae
const Foam::PtrList< Foam::mappedList< Foam::scalar> >&
Foam::conditionalMomentInversion::abscissae() const
{
    return abscissae_;
}


// ************************************************************************* //
