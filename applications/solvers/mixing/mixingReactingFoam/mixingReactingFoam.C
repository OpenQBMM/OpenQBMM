/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2026 Alberto Passalacqua
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
    mixingReactingFoam

Description
    Transient solver for the turbulent flow of a fluid whose composition and
    temperature are those of a mixing model: a mixture fraction, the
    compositions and the enthalpy conditioned on it, mixed and reacted in
    the environments of a quadrature, in the place of the transport
    equations of the species and of the energy.

    The solution follows rhoReactingBuoyantFoam: the continuity equation,
    the momentum equation and a pressure equation for p_rgh in a PIMPLE
    loop, on a thermo that holds its density. The thermo is selected at run
    time among those of that kind: a perfect gas for combustion, a liquid
    of constant density or a perfect fluid for the mixing of liquids, one
    fluid where nothing reacts and a reacting mixture where something does.

    The mixing model is solved once in a time step, right after the
    continuity equation and before the momentum equation of the first
    outer iteration. There the density obeys continuity with the mass flux
    that carries the moments, which is what keeps the moments of a
    mass-weighted distribution realizable and a distribution that
    integrates to one doing so. The model gives the thermo the mean
    composition and temperature of its environments, corrects it and sets
    its density to their mean density, which the pressure equation then
    reads: the flow responds to the heat the environments release. No
    energy equation is solved: with a model that gives the thermo no
    temperature, the temperature stays the one the case starts from.

    A mixing model of one coordinate, a mixture fraction alone, is
    transported by the volume flux, phi over the density of the faces, and
    is for a fluid of constant density.

    With frozenFlow on in the PIMPLE dictionary the flow is not solved: the
    velocity, the fluxes, the density and the turbulence stay what the case
    starts from, and only the mixing model advances, which is what
    mixingTransportFoam does.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "mixingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for the turbulent flow of a fluid whose composition"
        " and temperature are those of a mixing model, with reaction in the"
        " environments of its quadrature."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"

    turbulence->validate();

    #include "compressibleCourantNo.H"
    #include "realizableCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "realizableCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const bool frozenFlow =
            pimple.dict().getOrDefault<bool>("frozenFlow", false);

        if (!frozenFlow)
        {
            #include "rhoEqn.H"
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // The mixing model advances once in a time step, with the
            // density the continuity equation has just given
            if (pimple.firstIter())
            {
                phiByRho = phi/fvc::interpolate(rho);
                mixing->solve();
            }

            if (!frozenFlow)
            {
                #include "UEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }
        }

        if (!frozenFlow)
        {
            rho = thermo.rho();
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
