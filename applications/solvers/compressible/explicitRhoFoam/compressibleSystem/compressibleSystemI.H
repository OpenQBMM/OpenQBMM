/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 Alberto Passalacqua
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

const Foam::rhoThermo& Foam::compressibleSystem::thermo() const
{
    return thermoPtr_();
}


Foam::rhoThermo& Foam::compressibleSystem::thermo()
{
    return thermoPtr_();
}


const Foam::volScalarField& Foam::compressibleSystem::rho() const
{
    return rho_;
}


Foam::volScalarField& Foam::compressibleSystem::rho()
{
    return rho_;
}


const Foam::volScalarField& Foam::compressibleSystem::E() const
{
    return E_;
}


Foam::volScalarField& Foam::compressibleSystem::E()
{
    return E_;
}


const Foam::volScalarField& Foam::compressibleSystem::H() const
{
    return H_;
}


Foam::volScalarField& Foam::compressibleSystem::H()
{
    return H_;
}


const Foam::volVectorField& Foam::compressibleSystem::U() const
{
    return U_;
}


Foam::volVectorField& Foam::compressibleSystem::U()
{
    return U_;
}


const Foam::surfaceScalarField& Foam::compressibleSystem::phi() const
{
    return phiPtr_();
}


Foam::surfaceScalarField& Foam::compressibleSystem::phi()
{
    return phiPtr_();
}


const Foam::volScalarField& Foam::compressibleSystem::p() const
{
    return p_;
}


Foam::volScalarField& Foam::compressibleSystem::p()
{
    return p_;
}


const Foam::volVectorField& Foam::compressibleSystem::rhoU() const
{
    return rhoU_;
}


Foam::volVectorField& Foam::compressibleSystem::rhoU()
{
    return rhoU_;
}


const Foam::volScalarField& Foam::compressibleSystem::rhoE() const
{
    return rhoE_;
}


Foam::volScalarField& Foam::compressibleSystem::rhoE()
{
    return rhoE_;
}


const Foam::surfaceScalarField&
Foam::compressibleSystem::massFlux() const
{
    return massFlux_;
}


Foam::surfaceScalarField& Foam::compressibleSystem::massFlux()
{
    return massFlux_;
}


const Foam::surfaceVectorField&
Foam::compressibleSystem::momentumFlux() const
{
    return momentumFlux_;
}


Foam::surfaceVectorField& Foam::compressibleSystem::momentumFlux()
{
    return momentumFlux_;
}


const Foam::surfaceScalarField&
Foam::compressibleSystem::energyFlux() const
{
    return energyFlux_;
}


Foam::surfaceScalarField& Foam::compressibleSystem::energyFlux()
{
    return energyFlux_;
}


// ************************************************************************* //
