/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2019 Jeffrey Heylmun
    Copyright (C) 2020-2021 Alberto Passalacqua
-------------------------------------------------------------------------------
02-09-2019 Jeffrey Heylmun: Added Solution of moment trasporport as well as the
                            coupling of the continuous phase with the velocity
                            abscissae.
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

Application
    diluteVdfTransportFoam

Description
    Transient solver which couples a continuous phase with a dilute dispersed
    phase through force terms. Particle advection is done using quadrature based
    method of moments.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "constants.H"
#include "phaseModel.H"
#include "vdfPhaseModel.H"
#include "orderedPhasePair.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"
#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver which couples a continuous phase with a dilute "
        "dispersed phase through momentum exchange terms."
    );

    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "contErrs.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    Switch solveEnergy(phaseProperties.lookupOrDefault("solveEnergy", false));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phase2->solve();
        alpha1 = 1.0 - alpha2;
        alphaPhi1 = fvc::interpolate(alpha1)*phi1;
        alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
        K2 = 0.5*magSqr(U2);

        #include "momentumSources.H"

        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "correctContErrs.H"
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                phase1->correct();
                phase1->turbulence().correct();
            }
        }
        #include "vEqns.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
