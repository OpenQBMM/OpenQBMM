/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2016-06-30 Alberto Passalacqua: Implemented quadrature-based bubbly flow solver.
                                Implementation performed by J. C. Heylmun
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
    polydisperseBubbleFoam

Description
    Solver for polydisperse gas-liquid flows using a quadrature-based moment
    method to describe the gas phase.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvcSmooth.H"
#include "MULES.H"
#include "fvOptions.H"
#include "subCycle.H"

#include "phaseModel.H"
#include "pdPhaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"

#include "dragModel.H"
#include "liftModel.H"
#include "virtualMassModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"
#include "bubblePressureModel.H"
#include "BlendedInterfacialModel.H"
#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "createFields.H"

    #include "pU/createDDtU.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        {
            // Create interfacial forces and coefficients for use in velocity
            // moment transport as well as two fluid model
            #include "createInterfacialForces.H"

            // Transport moments with velocities relative to the mean gas
            // velocity
            phase1.relativeTransport();
            phi == fvc::interpolate(alpha1)*phi1 + fvc::interpolate(alpha2)*phi2;

            // Solve for mean phase velocities and gas volume fraction
            while (pimple.loop())
            {
                #include "alphas.H"

                #include "contErrs.H"
                #include "pU/DDtU.H"
                #include "updateInterfacialForces.H"

                #include "pU/UEqns.H"
                #include "pU/pEqn.H"
            }

            // Transport moments with mean gas velocity
            phase1.averageTransport(AEqns);

            phi ==
                fvc::interpolate(alpha1)*phi1
              + fvc::interpolate(alpha2)*phi2;

            if (pimple.turbCorr())
            {
                phase1.turbulence().correct();
                phase2.turbulence().correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
