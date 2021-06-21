/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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
    denseAGFoam

Description
    Solver for gas-particle flows with anisotropic Gaussian closures. It uses 
    a normal distribution of velocities in the dilute regime for more
    accurate transport. The 0th, 1st, and 2nd order velocity moments are used
    for transport. A switching function is also used to calculate the
    importance of dense and dilute solutions.

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
#include "MULES.H"
#include "subCycle.H"
#include "UniformField.H"

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
    argList::addNote
    (
        "Solver for gas-particle flows with anisotropic Gaussian closures."
    );
    
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

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        AGmodel.transportMoments();
        alpha2 = 1.0 - alpha1;
        alpha2.correctBoundaryConditions();

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            tmp<surfaceScalarField> h2fTmp(AGmodel.h2f());
            const surfaceScalarField& h2f = h2fTmp();

            #include "contErrs.H"

            #include "alphaEqn.H"
            fluid.correct();

			#include "pU/UEqns.H"
            if (solveP_rgh)
            {
                #include "pU/p_rghEqn.H"
            }
            else
            {
                #include "pU/pEqn.H"
            }

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
