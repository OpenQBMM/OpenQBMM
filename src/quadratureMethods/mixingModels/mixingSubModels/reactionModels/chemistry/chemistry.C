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

#include "chemistry.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "BasicChemistryModel.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixingSubModels
{
namespace reactionModels
{
    defineTypeNameAndDebug(chemistry, 0);

    addToRunTimeSelectionTable
    (
        reactionModel,
        chemistry,
        dictionary
    );
}
}
}

const Foam::Enum
<
    Foam::mixingSubModels::reactionModels::chemistry::temperatureClosure
>
Foam::mixingSubModels::reactionModels::chemistry::temperatureClosureNames
({
    { temperatureClosure::isothermal, "isothermal" },
    { temperatureClosure::linearEnthalpy, "linearEnthalpy" },
    { temperatureClosure::conditionalEnthalpy, "conditionalEnthalpy" },
});


// * * * * * * * * * * * * * * isothermalSystem  * * * * * * * * * * * * * * //

void Foam::mixingSubModels::reactionModels::chemistry::isothermalSystem
::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    system_.derivatives(x, y, dydx);
    dydx[TIndex_] = 0;
}


void Foam::mixingSubModels::reactionModels::chemistry::isothermalSystem
::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    system_.jacobian(x, y, dfdx, dfdy);
    dfdx[TIndex_] = 0;

    for (label j = 0; j < dfdy.n(); j++)
    {
        dfdy(TIndex_, j) = 0;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::basicSpecieMixture&
Foam::mixingSubModels::reactionModels::chemistry::compositionOf
(
    basicThermo& thermo
)
{
    if (isA<psiReactionThermo>(thermo))
    {
        return refCast<psiReactionThermo>(thermo).composition();
    }
    else if (isA<rhoReactionThermo>(thermo))
    {
        return refCast<rhoReactionThermo>(thermo).composition();
    }

    FatalErrorInFunction
        << "The chemistry reaction model needs a thermo with species and "
        << "reactions, a psiReactionThermo or a rhoReactionThermo with a "
        << "reactingMixture; the thermo of the case is "
        << thermo.type() << nl
        << exit(FatalError);

    return refCast<psiReactionThermo>(thermo).composition();
}


Foam::autoPtr<Foam::basicChemistryModel>
Foam::mixingSubModels::reactionModels::chemistry::newChemistry
(
    basicThermo& thermo
)
{
    if (isA<psiReactionThermo>(thermo))
    {
        return autoPtr<basicChemistryModel>
        (
            basicChemistryModel::New<BasicChemistryModel<psiReactionThermo>>
            (
                refCast<psiReactionThermo>(thermo)
            ).ptr()
        );
    }

    return autoPtr<basicChemistryModel>
    (
        basicChemistryModel::New<BasicChemistryModel<rhoReactionThermo>>
        (
            refCast<rhoReactionThermo>(thermo)
        ).ptr()
    );
}


const Foam::ODESystem&
Foam::mixingSubModels::reactionModels::chemistry::systemOf
(
    const basicChemistryModel& model
)
{
    const ODESystem* systemPtr = dynamic_cast<const ODESystem*>(&model);

    if (!systemPtr)
    {
        FatalErrorInFunction
            << "The chemistry model " << model.type()
            << " does not give the derivatives of the state of an "
            << "environment, which the chemistry reaction model integrates."
            << nl << exit(FatalError);
    }

    return *systemPtr;
}


void Foam::mixingSubModels::reactionModels::chemistry::readStream
(
    const dictionary& dict,
    scalar& T,
    scalarField& Y
) const
{
    T = dict.get<scalar>("T");

    const dictionary& YDict = dict.subDict("Y");
    const speciesTable& species = composition_.species();

    Y = scalarField(species.size(), Zero);

    for (const entry& e : YDict)
    {
        const label speciei = species.find(e.keyword());

        if (speciei < 0)
        {
            FatalIOErrorInFunction(YDict)
                << e.keyword() << " is not a species of the thermo, whose "
                << "species are " << species << nl
                << exit(FatalIOError);
        }

        Y[speciei] = YDict.get<scalar>(e.keyword());
    }

    const scalar YSum = sum(Y);

    if (mag(YSum - 1) > 1e-6)
    {
        FatalIOErrorInFunction(YDict)
            << "The mass fractions of the stream sum to " << YSum
            << ", not to one." << nl
            << exit(FatalIOError);
    }

    Y /= YSum;
}


Foam::scalar Foam::mixingSubModels::reactionModels::chemistry::enthalpy
(
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    scalar h = 0;

    forAll(Y, speciei)
    {
        h += Y[speciei]*composition_.Ha(speciei, p, T);
    }

    return h;
}


Foam::scalar Foam::mixingSubModels::reactionModels::chemistry::Cp
(
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    scalar Cp = 0;

    forAll(Y, speciei)
    {
        Cp += Y[speciei]*composition_.Cp(speciei, p, T);
    }

    return Cp;
}


Foam::scalar Foam::mixingSubModels::reactionModels::chemistry::density
(
    const scalarField& Y,
    const scalar p,
    const scalar T
) const
{
    // The species share the volume: the specific volume of the mixture is
    // the mass-weighted sum of theirs
    scalar v = 0;

    forAll(Y, speciei)
    {
        if (Y[speciei] > 0)
        {
            v += Y[speciei]/composition_.rho(speciei, p, T);
        }
    }

    return 1/v;
}


Foam::scalar Foam::mixingSubModels::reactionModels::chemistry::nodeTemperature
(
    const scalar xi,
    const scalarField& Y,
    const scalar hNode,
    const scalar Tcell,
    const scalar p
) const
{
    if (temperature_ == temperatureClosure::isothermal)
    {
        return Tcell;
    }

    // The enthalpy of the environment: the one its node carries, or the
    // one linear in the mixture fraction between the streams
    scalar h = hNode;

    if (temperature_ == temperatureClosure::linearEnthalpy)
    {
        const scalar h0 = enthalpy(Y0_, p, T0_);
        const scalar h1 = enthalpy(Y1_, p, T1_);
        h = h0 + xi*(h1 - h0);
    }

    // The temperature at which the composition has that enthalpy, by
    // Newton iteration from the temperature of the cell
    scalar T = Tcell;

    for (label iter = 0; iter < 100; iter++)
    {
        const scalar dT = (enthalpy(Y, p, T) - h)/Cp(Y, p, T);
        T = max(T - dT, SMALL);

        if (mag(dT) < 1e-10*T)
        {
            return T;
        }
    }

    FatalErrorInFunction
        << "The temperature of an environment of mixture fraction " << xi
        << " and composition " << Y << " with enthalpy " << h
        << " did not converge: " << T << " K after 100 iterations." << nl
        << exit(FatalError);

    return T;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingSubModels::reactionModels::chemistry::chemistry
(
    const dictionary& dict,
    const velocityQuadratureApproximation& quadrature,
    const label nCompositions
)
:
    reactionModel(dict, quadrature, nCompositions),
    thermo_(mesh_.lookupObjectRef<basicThermo>(basicThermo::dictName)),
    composition_(compositionOf(thermo_)),
    chemistry_(newChemistry(thermo_)),
    system_(systemOf(chemistry_())),
    temperature_(temperatureClosureNames.get("temperature", dict)),
    isothermalSystem_
    (
        temperature_ == temperatureClosure::isothermal
      ? new isothermalSystem(system_, system_.nEqns() - 2)
      : nullptr
    ),
    odeSolver_
    (
        ODESolver::New
        (
            temperature_ == temperatureClosure::isothermal
          ? static_cast<const ODESystem&>(isothermalSystem_())
          : system_,
            chemistry_->subDict("odeCoeffs")
        )
    ),
    T0_(0),
    Y0_(),
    T1_(0),
    Y1_(),
    Treact_(chemistry_->getOrDefault<scalar>("Treact", 0)),
    deltaTChemIni_(chemistry_->get<scalar>("initialChemicalTimeStep")),
    deltaTChemMax_
    (
        chemistry_->getOrDefault<scalar>("maxChemicalTimeStep", GREAT)
    ),
    subDeltaT_
    (
        quadrature.nodes().size(),
        scalarField(mesh_.nCells(), deltaTChemIni_)
    ),
    T_(quadrature.nodes().size()),
    cTp_(system_.nEqns(), Zero),
    Y_(nCompositions, Zero)
{
    const label nSpecie = composition_.species().size();

    forAll(T_, nodei)
    {
        T_.set
        (
            nodei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("T", quadrature.nodes()[nodei].name()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo_.T(),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    if (nSpecie != nCompositions)
    {
        FatalIOErrorInFunction(dict)
            << "The quadrature carries " << nCompositions
            << " compositions conditioned on the mixture fraction, and the "
            << "thermo has " << nSpecie << " species: "
            << composition_.species() << nl
            << "    The compositions are the species, in the order of the "
            << "species list." << nl
            << exit(FatalIOError);
    }

    if (system_.nEqns() != nSpecie + 2)
    {
        FatalErrorInFunction
            << "The chemistry model " << chemistry_->type()
            << " integrates " << system_.nEqns() << " equations, not the "
            << nSpecie << " species, the temperature and the pressure."
            << nl << exit(FatalError);
    }

    if (temperature_ == temperatureClosure::linearEnthalpy)
    {
        const dictionary& streams = dict.subDict("streams");

        readStream(streams.subDict("xi0"), T0_, Y0_);
        readStream(streams.subDict("xi1"), T1_, Y1_);
    }

    Info<< "    " << nSpecie << " species react in every environment, "
        << "temperature " << temperatureClosureNames[temperature_] << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingSubModels::reactionModels::chemistry::~chemistry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mixingSubModels::reactionModels::chemistry::reacts() const
{
    return chemistry_->chemistry();
}


void Foam::mixingSubModels::reactionModels::chemistry::react
(
    velocityQuadratureApproximation& quadrature
)
{
    if (!reacts())
    {
        return;
    }

    const label nSpecie = nCompositions_;
    const scalar deltaT = mesh_.time().deltaTValue();
    const scalarField& p = thermo_.p().primitiveField();
    const scalarField& T = thermo_.T().primitiveField();

    mappedPtrList<volVelocityNode>& nodes = quadrature.nodes();
    const label nNodes = nodes.size();

    // The abscissae of a node hold the mixture fraction first, then the
    // compositions, and last the enthalpy when the distribution carries it
    const bool enthalpyCarried =
        temperature_ == temperatureClosure::conditionalEnthalpy;

    if
    (
        nodes[0].scalarIndexes().size()
     != 1 + nSpecie + (enthalpyCarried ? 1 : 0)
    )
    {
        FatalErrorInFunction
            << "The nodes of the quadrature carry "
            << nodes[0].scalarIndexes().size()
            << " scalars, not the mixture fraction and " << nSpecie
            << " compositions"
            << (enthalpyCarried ? " and the enthalpy." : ".") << nl
            << exit(FatalError);
    }

    if (subDeltaT_.size() != nNodes)
    {
        FatalErrorInFunction
            << "The quadrature has " << nNodes << " nodes, and the model "
            << "was built for " << subDeltaT_.size() << "." << nl
            << exit(FatalError);
    }

    scalarField meanT(T);
    scalarList meanY(nSpecie, Zero);

    forAll(T, celli)
    {
        const scalar pi = p[celli];

        scalar wSum = 0;
        scalar wT = 0;
        meanY = Zero;

        forAll(nodes, nodei)
        {
            volVelocityNode& node = nodes[nodei];
            const scalar w = node.weight()[celli];

            // A node carrying nothing carries no composition to react
            if (w <= SMALL)
            {
                T_[nodei][celli] = T[celli];
                continue;
            }

            const scalar xi = node.abscissae()[0][celli];

            // The composition of the environment, within its bounds and
            // summing to one
            scalar YSum = 0;

            for (label j = 0; j < nSpecie; j++)
            {
                Y_[j] = min(max(node.abscissae()[1 + j][celli], 0), 1);
                YSum += Y_[j];
            }

            if (YSum <= 0)
            {
                T_[nodei][celli] = T[celli];
                continue;
            }

            Y_ /= YSum;

            const scalar hi =
                enthalpyCarried
              ? node.abscissae()[1 + nSpecie][celli]
              : 0;

            scalar Ti = nodeTemperature(xi, Y_, hi, T[celli], pi);

            if (Ti > Treact_)
            {
                const scalar rhoi = density(Y_, pi, Ti);

                for (label j = 0; j < nSpecie; j++)
                {
                    cTp_[j] = rhoi*Y_[j]/composition_.W(j);
                }

                cTp_[nSpecie] = Ti;
                cTp_[nSpecie + 1] = pi;

                scalar& subDeltaT = subDeltaT_[nodei][celli];

                if (debug)
                {
                    Info<< "chemistry: cell " << celli << " node " << nodei
                        << " w " << w << " xi " << xi << " Y " << Y_
                        << " T " << Ti << " rho " << rhoi
                        << " c " << cTp_ << " subDeltaT " << subDeltaT
                        << endl;
                }

                odeSolver_->solve(0, deltaT, cTp_, subDeltaT);

                if (debug)
                {
                    Info<< "chemistry: reacted to " << cTp_ << endl;
                }
                subDeltaT = min(subDeltaT, deltaTChemMax_);

                // Back to mass fractions, normalised by the mass the
                // concentrations carry
                scalar mass = 0;

                for (label j = 0; j < nSpecie; j++)
                {
                    cTp_[j] = max(cTp_[j], 0);
                    mass += cTp_[j]*composition_.W(j);
                }

                for (label j = 0; j < nSpecie; j++)
                {
                    Y_[j] = cTp_[j]*composition_.W(j)/mass;
                    node.abscissae()[1 + j][celli] = Y_[j];
                }

                Ti = cTp_[nSpecie];
            }

            T_[nodei][celli] = Ti;

            wSum += w;
            wT += w*Ti;

            for (label j = 0; j < nSpecie; j++)
            {
                meanY[j] += w*Y_[j];
            }
        }

        quadrature.updateLocalMoments(celli);

        if (wSum > 0)
        {
            for (label j = 0; j < nSpecie; j++)
            {
                composition_.Y(j)[celli] = meanY[j]/wSum;
            }

            meanT[celli] = wT/wSum;
        }
    }

    forAll(quadrature.moments(), mi)
    {
        quadrature.moments()[mi].correctBoundaryConditions();
    }

    quadrature.updateBoundaryQuadrature();

    for (label j = 0; j < nSpecie; j++)
    {
        composition_.Y(j).correctBoundaryConditions();
    }

    forAll(T_, nodei)
    {
        T_[nodei].correctBoundaryConditions();
    }

    if (temperature_ != temperatureClosure::isothermal)
    {
        thermo_.he().primitiveFieldRef() =
            thermo_.he(p, meanT, identity(mesh_.nCells()));
    }
}


// ************************************************************************* //
