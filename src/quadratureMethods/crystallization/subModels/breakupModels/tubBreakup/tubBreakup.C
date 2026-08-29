/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
#include "tubBreakup.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulenceModel.H"

namespace Foam
{
namespace populationBalanceSubModels
{
namespace breakupKernels
{
    defineTypeNameAndDebug(tubBreakup, 0);
    addToRunTimeSelectionTable
    (
        breakupKernel,
        tubBreakup,
        dictionary
    );
}
}
}

Foam::populationBalanceSubModels::breakupKernels::tubBreakup::tubBreakup
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    breakupKernel(dict, mesh),
    continuousPhase_(dict.getOrDefault<word>("continuousPhase", word::null)),
    des_("des", dimDensity, dict),
    vis_("vis", dimDynamicViscosity, dict),
    Df_("Df", dimless, dict),
    R0_("R0", dimLength, dict),
    useConstantEpsilon_(dict.getOrDefault<bool>("useConstantEpsilon", false)),
    constantEpsilon_("constantEpsilon", dimVelocity*dimVelocity/dimTime, dict),
    flTurb_(nullptr),
    epsilon_(),
    gField_
    (
        IOobject
        (
            "breakupRate",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    )
{
    if (!useConstantEpsilon_)
    {
        const word turbName = IOobject::groupName
        (
            turbulenceModel::propertiesName,
            continuousPhase_
        );

        if (mesh_.foundObject<turbulenceModel>(turbName))
        {
            flTurb_ = &mesh_.lookupObject<turbulenceModel>(turbName);
            epsilon_ = flTurb_->epsilon();
        }
        else
        {
            FatalErrorInFunction
                << "No valid turbulence model found for phase "
                << continuousPhase_ << nl
                << exit(FatalError);
        }
    }
}

Foam::populationBalanceSubModels::breakupKernels::tubBreakup::~tubBreakup()
{}

Foam::scalar
Foam::populationBalanceSubModels::breakupKernels::tubBreakup::Kb
(
    const scalar& abscissa,
    const label celli,
    const label environment
) const
{
    const scalar nu = vis_.value() / des_.value();
    //const scalar R = pow(2, -3.0/Df_.value()) * pow(R0_.value(), 1.0 - 3.0/Df_.value()) 
               //* pow(abscissa, 3.0/Df_.value());
    
    const scalar localEpsilon = useConstantEpsilon_ ? 
        constantEpsilon_.value() : 
        epsilon_()[celli];
    
    const scalar G = sqrt(localEpsilon / nu);
    const scalar b = 0.2159 * pow(G, 1.85);
    const scalar g = b * pow(abscissa, 2);
    
    // 更新 gField_
    const_cast<volScalarField&>(gField_)[celli] = g;

    return g;
}

// ************************************************************************* //
