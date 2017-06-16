/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipUgation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be usefUg, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICUgAR PURPOSE.  See the GNU General Public License
    for more details.

    You shoUgd have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", 0.01);

scalar	CoNumParticleKinetics = AGmodel.maxUxDx()*runTime.deltaTValue();

scalarField sumPhi =  fvc::surfaceSum(mag(phi1))().internalField();
scalar	CoNumParticlePhi =
    0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

scalar CoNumParticle = max(CoNumParticleKinetics,CoNumParticlePhi);

sumPhi =  fvc::surfaceSum(mag(phi2))().internalField();
scalar	CoNumGas = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

scalar	CoNum = max(CoNumParticle,CoNumGas);

scalar maxDeltaTFact = maxCo/(CoNum + SMALL);

scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

runTime.setDeltaT
(
    min
    (
        deltaTFact*runTime.deltaTValue(),
        maxDeltaT
    )
);
Info << "CoNumPartKin  = " <<  CoNumParticleKinetics  << endl;
Info << "CoNumPartPhi  = " <<  CoNumParticlePhi  << endl;
Info << "CoNumGas      = " <<  CoNumGas  << endl;
Info << "deltaT = " <<  runTime.deltaTValue() << nl << endl;

// ************************************************************************* //
