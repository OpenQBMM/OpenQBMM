/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    twoPhaseEulerFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "kineticTheoryModel.H"
#include "fixedValueFvsPatchFields.H"

#include "fixedFluxPressureFvPatchScalarField.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
void explicitSolve
(
    volScalarField& psi,
    const surfaceScalarField& phiPsi,
    const scalar& deltaT
)
{
    scalarField& psiIf = psi;
    const scalarField& psi0 = psi.oldTime();

    psiIf = 0.0;
    fvc::surfaceIntegrate(psiIf, phiPsi);

    psiIf = psi0 - psiIf*deltaT;

    psi.correctBoundaryConditions();

}
*/

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    Switch implicitPhasePressure
    (
        mesh.solverDict(alpha1.name()).lookupOrDefault<Switch>
        (
            "implicitPhasePressure", false
        )
    );

    #include "pU/createDDtU.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"
        if (adjustTimeStep)
        {
            runTime.setDeltaT
            (
                min
                (
                    runTime.deltaT(),
                    AGmodel.maxUxDx()*runTime.deltaT()
                )
            );
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        AGmodel.transportMoments();

        //  Update alpha2 and set ddt(alpha1)
        alpha2 = 1.0 - alpha1;
        volScalarField ddtAlpha1Dilute(fvc::ddt(alpha1));
        volVectorField ddtU1Dilute(fvc::ddt(alpha1, rho1, U1));


        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "contErrs.H"
            {
                //  Store old phi1 so that the dense flux is used to compute
                //  the volume fraction transport
                surfaceScalarField phiOld = phi1;

                surfaceScalarField pPrimeByA = fluid.pPrimeByA()();
                phi1 = AGmodel.hydrodynamicScalef
                (
                    phi1
                  + pPrimeByA*fvc::snGrad(alpha1, "bounded")*mesh.magSf()
                );

                alphaPhi1 = fvc::interpolate(alpha1)*phi1;
                alphaRhoPhi1 = phase1.alphaPhi()*fvc::interpolate(rho1);

                fluid.pPrimeByA().ref() =
                    AGmodel.hydrodynamicScalef(pPrimeByA);

                fluid.solve();
                fluid.correct();

                phi1 = phiOld;
            }

			#include "pU/UEqns.H"
            #include "pU/pEqn.H"
            #include "pU/DDtU.H"

            if (pimple.turbCorr())
            {
				fluid.correctTurbulence();
            }
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime.write();

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //