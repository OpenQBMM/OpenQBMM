/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

const Foam::PtrList<Foam::volVectorField>&
Foam::monoKineticQuadratureApproximation::velocities() const
{
    return velocityAbscissae_;
}

Foam::PtrList<Foam::volVectorField>&
Foam::monoKineticQuadratureApproximation::velocities()
{
    return velocityAbscissae_;
}

const Foam::PtrList<Foam::volVectorField>&
Foam::monoKineticQuadratureApproximation::velocityMoments() const
{
    return velocityMoments_;
}

Foam::PtrList<Foam::volVectorField>&
Foam::monoKineticQuadratureApproximation::velocityMoments()
{
    return velocityMoments_;
}

const Foam::PtrList<Foam::surfaceNode>&
Foam::monoKineticQuadratureApproximation::nodesOwn() const
{
    return nodesOwn_();
}

const Foam::PtrList<Foam::surfaceVectorField>&
Foam::monoKineticQuadratureApproximation::velocitiesOwn() const
{
    return velocitiesOwn_;
}

const Foam::PtrList<Foam::surfaceNode>&
Foam::monoKineticQuadratureApproximation::nodesNei() const
{
    return nodesNei_;
}

const Foam::PtrList<Foam::surfaceVectorField>&
Foam::monoKineticQuadratureApproximation::velocitiesNei() const
{
    return velocitiesNei_;
}

// ************************************************************************* //
