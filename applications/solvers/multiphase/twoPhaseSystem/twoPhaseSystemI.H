/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

const Foam::fvMesh& Foam::twoPhaseSystem::mesh() const
{
    return mesh_;
}


const Foam::phaseModel& Foam::twoPhaseSystem::phase1() const
{
    return phase1_();
}


Foam::phaseModel& Foam::twoPhaseSystem::phase1()
{
    return phase1_();
}


const Foam::phaseModel& Foam::twoPhaseSystem::phase2() const
{
    return phase2_();
}


Foam::phaseModel& Foam::twoPhaseSystem::phase2()
{
    return phase2_();
}


const Foam::phaseModel& Foam::twoPhaseSystem::otherPhase
(
    const phaseModel& phase
) const
{
    if (&phase == &phase1_())
    {
        return phase2_();
    }
    else
    {
        return phase1_();
    }
}

const Foam::phasePair& Foam::twoPhaseSystem::pair() const
{
    return pair_();
}


const Foam::orderedPhasePair& Foam::twoPhaseSystem::pair1In2() const
{
    return pair1In2_();
}


const Foam::orderedPhasePair& Foam::twoPhaseSystem::pair2In1() const
{
    return pair2In1_();
}


const Foam::surfaceScalarField& Foam::twoPhaseSystem::phi() const
{
    return phi_;
}


Foam::surfaceScalarField& Foam::twoPhaseSystem::phi()
{
    return phi_;
}


const Foam::volScalarField& Foam::twoPhaseSystem::dgdt() const
{
    return dgdt_;
}


Foam::volScalarField& Foam::twoPhaseSystem::dgdt()
{
    return dgdt_;
}


Foam::tmp<Foam::surfaceScalarField>& Foam::twoPhaseSystem::pPrimeByA()
{
    return pPrimeByA_;
}

// ************************************************************************* //
