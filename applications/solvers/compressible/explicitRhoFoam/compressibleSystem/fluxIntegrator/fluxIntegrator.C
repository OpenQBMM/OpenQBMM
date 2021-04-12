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

#include "fluxIntegrator.H"


// * * * * * * * * * * * * * * Private Data Members * * * * * * * * * * * * * //

void Foam::fluxIntegrator::setCoeffs
(
    boolList& storeFields,
    boolList& storeDeltas
)
{
    List<scalarList> c = butcherTable_->conservedVariablesCoeffs();
    List<scalarList> f = butcherTable_->fluxCoeffs();

    for (label i = 0; i < c.size(); i++)
    {
        for (label j = 0; j < c[i].size() - 1; j++)
        {
            if (mag(c[i][j]) > SMALL)
            {
                storeFields[j] = true;
            }

            if (mag(f[i][j]) > SMALL)
            {
                storeDeltas[j] = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxIntegrator::fluxIntegrator(compressibleSystem& fluid)
:
    fluid_(fluid),
    butcherTable_(ButcherTable::New(fluid.rho().mesh()))
{
    boolList storeFields(butcherTable_->nSteps(), false);
    boolList storeDeltas(butcherTable_->nSteps(), false);
    setCoeffs(storeFields, storeDeltas);
    fluid_.setNSteps(storeFields, storeDeltas);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxIntegrator::~fluxIntegrator()
{}


// * * * * * * * * * * * * * *  Public Functinos * * * * * * * * * * * * * * //

void Foam::fluxIntegrator::integrateFluxes(const dimensionedVector& g)
{
    List<scalarList> c = butcherTable_->conservedVariablesCoeffs();
    List<scalarList> f = butcherTable_->fluxCoeffs();

    const dimensionedScalar& deltaT = fluid_.rho().time().deltaT();
    
    fluid_.calcConservativeVariables();
    
    for (label stepi = 0; stepi < butcherTable_->nSteps(); stepi++)
    {
        fluid_.updateFluxes();

        fluid_.advect
        (
            stepi,
            c[stepi],
            f[stepi],
            deltaT,
            g
        );

        fluid_.calcPrimitiveVariables();
    }
}

// ************************************************************************* //
