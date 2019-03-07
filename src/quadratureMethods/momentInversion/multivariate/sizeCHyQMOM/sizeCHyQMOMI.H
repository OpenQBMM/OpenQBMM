/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the total number of moments
Foam::label Foam::sizeCHyQMOM::nMoments()
{
    return nMoments_;
}

//- Return the list of nNodes
Foam::label Foam::sizeCHyQMOM::nNodes()
{
    return nNodes_;
}

//- Return mapped list of weights
const Foam::mappedList<Foam::scalar>& Foam::sizeCHyQMOM::weights() const
{
    return weights_;
}

//- Return mapped list of abscissae
const Foam::mappedList<Foam::scalarList>& Foam::sizeCHyQMOM::abscissae() const
{
    return abscissae_;
}


// ************************************************************************* //
