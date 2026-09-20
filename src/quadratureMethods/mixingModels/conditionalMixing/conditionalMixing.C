/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "conditionalMixing.H"
#include "basicThermo.H"
#include "basicSpecieMixture.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace mixingModels
{
    defineTypeNameAndDebug(conditionalMixing, 0);

    addToRunTimeSelectionTable
    (
        mixingModel,
        conditionalMixing,
        dictionary
    );
}
}

//- The moment fields of a set, as the mixing kernels read them
static UPtrList<const volScalarField> momentFieldsOf
(
    const volVelocityMomentFieldSet& moments
)
{
    UPtrList<const volScalarField> fields(moments.size());

    forAll(fields, momenti)
    {
        fields.set(momenti, &moments[momenti]);
    }

    return fields;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::label
Foam::PDFTransportModels::mixingModels::conditionalMixing::nCompositions
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    label nCompositions = -1;

    if (dict.readIfPresent("nCompositions", nCompositions))
    {
        if (nCompositions < 0)
        {
            FatalIOErrorInFunction(dict)
                << "The number of compositions conditioned on the mixture "
                << "fraction cannot be negative." << nl
                << exit(FatalIOError);
        }

        return nCompositions;
    }

    // Without nCompositions the compositions are the species of the thermo
    const basicThermo* thermoPtr =
        mesh.cfindObject<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        if (isA<psiReactionThermo>(*thermoPtr))
        {
            return
                refCast<const psiReactionThermo>(*thermoPtr)
               .composition().species().size();
        }
        else if (isA<rhoReactionThermo>(*thermoPtr))
        {
            return
                refCast<const rhoReactionThermo>(*thermoPtr)
               .composition().species().size();
        }
    }

    FatalIOErrorInFunction(dict)
        << "The number of compositions conditioned on the mixture fraction "
        << "is nCompositions, or the number of species of a thermo with "
        << "species, which the case does not have." << nl
        << exit(FatalIOError);

    return nCompositions;
}


Foam::List<Foam::supportType>
Foam::PDFTransportModels::mixingModels::conditionalMixing::supports
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    // The mixture fraction and every composition are mass fractions
    return
        List<supportType>
        (
            1 + nCompositions(dict, mesh),
            supportType::ZeroOne
        );
}


Foam::PDFTransportModels::mixingModels::conditionalMixing::conditionalMixing
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    velocityPDFTransportModel
    (
        name,
        dict,
        phi.mesh(),
        supports(dict, phi.mesh())
    ),
    mixingModel(name, dict, phi),
    odeType(phi.mesh(), dict),
    name_(name),
    momentFields_(momentFieldsOf((*this).quadrature().moments())),
    mixingKernel_
    (
        Foam::mixingSubModels::mixingKernel::New
        (
            dict.subDict("mixingKernel"),
            phi.mesh(),
            (*this).quadrature().momentOrders(),
            momentFields_
        )
    ),
    diffusionModel_
    (
        Foam::mixingSubModels::mixingDiffusionModel::New
        (
            dict.subDict("diffusionModel")
        )
    ),
    reactionModel_
    (
        Foam::mixingSubModels::reactionModel::New
        (
            dict.subOrEmptyDict("reaction"),
            quadrature_,
            nCompositions(dict, phi.mesh())
        )
    )
{
    const label nDimensions = quadrature().momentOrders()[0].size();
    const label nCompositions = reactionModel_->nCompositions();

    if (nDimensions != 1 + nCompositions)
    {
        FatalIOErrorInFunction(dict)
            << "The quadrature of " << name_ << " carries " << nDimensions
            << " dimensions, and the model mixes a mixture fraction with "
            << nCompositions << " compositions." << nl
            << "    Give nScalars of the inversion the number of "
            << "compositions." << nl
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::mixingModels::conditionalMixing::~conditionalMixing()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PDFTransportModels::mixingModels::conditionalMixing
::updateImplicitMomentSource()
{}


Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::mixingModels::conditionalMixing
::implicitMomentSource(const volVelocityMoment& moment)
{
    return diffusionModel_->momentDiff(moment);
}


void Foam::PDFTransportModels::mixingModels::conditionalMixing
::explicitMomentSource()
{
    odeType::solve(quadrature_, 0);
}


bool Foam::PDFTransportModels::mixingModels::conditionalMixing
::solveMomentSources() const
{
    return odeType::solveSources_;
}


bool Foam::PDFTransportModels::mixingModels::conditionalMixing
::solveMomentOde() const
{
    return odeType::solveOde_;
}


void Foam::PDFTransportModels::mixingModels::conditionalMixing
::updateCellMomentSource(const label)
{}


Foam::scalar
Foam::PDFTransportModels::mixingModels::conditionalMixing::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation& quadrature,
    const label environment
)
{
    return mixingKernel_->mixingSource(momentOrder, celli, environment);
}


Foam::scalar
Foam::PDFTransportModels::mixingModels::conditionalMixing::realizableCo()
{
    return velocityPDFTransportModel::realizableCo();
}


void Foam::PDFTransportModels::mixingModels::conditionalMixing::solve()
{
    velocityPDFTransportModel::solve();

    if (reactionModel_->reacts())
    {
        reactionModel_->react(quadrature_);
    }
}

// ************************************************************************* //
