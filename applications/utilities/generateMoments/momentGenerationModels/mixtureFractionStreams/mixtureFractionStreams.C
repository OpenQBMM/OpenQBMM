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

#include "mixtureFractionStreams.H"
#include "IOdictionary.H"
#include "IFstream.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{
    defineTypeNameAndDebug(mixtureFractionStreams, 0);

    addToRunTimeSelectionTable
    (
        momentGenerationModel,
        mixtureFractionStreams,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList
Foam::momentGenerationSubModels::mixtureFractionStreams::readSpecies
(
    const dictionary& dict
) const
{
    if (dict.found("species"))
    {
        return dict.get<wordList>("species");
    }

    const IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    if (!thermoDict.found("foamChemistryFile"))
    {
        FatalIOErrorInFunction(dict)
            << "The compositions are the species of the case, which are "
            << "read from the foamChemistryFile of thermophysicalProperties "
            << "when the dictionary has no species entry; the case names no "
            << "such file." << nl
            << exit(FatalIOError);
    }

    fileName chemistryFile(thermoDict.get<fileName>("foamChemistryFile"));
    chemistryFile.expand();

    IFstream chemistryStream(chemistryFile);

    if (!chemistryStream.good())
    {
        FatalIOErrorInFunction(thermoDict)
            << "Cannot read the species from " << chemistryFile << nl
            << exit(FatalIOError);
    }

    return dictionary(chemistryStream).get<wordList>("species");
}


Foam::scalarField
Foam::momentGenerationSubModels::mixtureFractionStreams::readStream
(
    const dictionary& dict
) const
{
    const dictionary& YDict = dict.subDict("Y");

    scalarField Y(species_.size(), Zero);

    for (const entry& e : YDict)
    {
        const label speciei = species_.find(e.keyword());

        if (speciei < 0)
        {
            FatalIOErrorInFunction(YDict)
                << e.keyword() << " is not a species of the case, whose "
                << "species are " << species_ << nl
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

    return Y;
}


void
Foam::momentGenerationSubModels::mixtureFractionStreams::readStreamEnthalpies
(
    const dictionary& streams
)
{
    // The thermo reads the mass fraction of every species, which must sum
    // to one: where the case does not have them yet, they are written with
    // the composition of the first stream, to be replaced by the means
    forAll(species_, j)
    {
        IOobject header
        (
            species_[j],
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        if (!header.typeHeaderOk<volScalarField>(true))
        {
            volScalarField
            (
                header,
                mesh_,
                dimensionedScalar(species_[j], dimless, Y0_[j])
            ).write();
        }
    }

    thermoPtr_ = fluidThermo::New(mesh_);

    const basicSpecieMixture* compositionPtr = nullptr;

    if (isA<psiReactionThermo>(thermoPtr_()))
    {
        compositionPtr =
            &refCast<psiReactionThermo>(thermoPtr_()).composition();
    }
    else if (isA<rhoReactionThermo>(thermoPtr_()))
    {
        compositionPtr =
            &refCast<rhoReactionThermo>(thermoPtr_()).composition();
    }
    else
    {
        FatalErrorInFunction
            << "The enthalpies of the streams need a thermo with species, "
            << "and the thermo of the case is " << thermoPtr_->type() << nl
            << exit(FatalError);
    }

    const basicSpecieMixture& composition = *compositionPtr;

    if (composition.species() != species_)
    {
        FatalErrorInFunction
            << "The species of the thermo, " << composition.species()
            << ", are not the species of the moments, " << species_ << nl
            << exit(FatalError);
    }

    const scalar p = gAverage(thermoPtr_->p().primitiveField());
    const scalar T0 = streams.subDict("xi0").get<scalar>("T");
    const scalar T1 = streams.subDict("xi1").get<scalar>("T");

    h0_ = 0;
    h1_ = 0;

    forAll(species_, j)
    {
        h0_ += Y0_[j]*composition.Ha(j, p, T0);
        h1_ += Y1_[j]*composition.Ha(j, p, T1);
    }

    Info<< "Absolute enthalpy of the streams: " << h0_ << " J/kg at xi = 0, "
        << h1_ << " J/kg at xi = 1" << endl;
}


Foam::scalarList
Foam::momentGenerationSubModels::mixtureFractionStreams
::mixtureFractionMoments
(
    const dictionary& dict,
    const label maxOrder
) const
{
    scalarList m(maxOrder + 1, Zero);
    m[0] = 1;

    const word distribution(dict.get<word>("distribution"));

    if (distribution == "beta")
    {
        const scalar mean = dict.get<scalar>("mean");
        const scalar variance = dict.get<scalar>("variance");

        if
        (
            mean <= 0 || mean >= 1
         || variance <= 0 || variance >= mean*(1 - mean)
        )
        {
            FatalIOErrorInFunction(dict)
                << "A beta distribution has 0 < mean < 1 and 0 < variance "
                << "< mean (1 - mean) = " << mean*(1 - mean)
                << "; the mean is " << mean << " and the variance "
                << variance << "." << nl
                << exit(FatalIOError);
        }

        const scalar c = mean*(1 - mean)/variance - 1;
        const scalar a = mean*c;
        const scalar b = (1 - mean)*c;

        for (label k = 1; k <= maxOrder; k++)
        {
            m[k] = m[k - 1]*(a + k - 1)/(a + b + k - 1);
        }
    }
    else if (distribution == "delta")
    {
        const scalar value = dict.get<scalar>("value");

        if (value < 0 || value > 1)
        {
            FatalIOErrorInFunction(dict)
                << "A mixture fraction is within [0, 1], and the value is "
                << value << "." << nl
                << exit(FatalIOError);
        }

        for (label k = 1; k <= maxOrder; k++)
        {
            m[k] = m[k - 1]*value;
        }
    }
    else if (distribution == "nodes")
    {
        const scalarList weights(dict.get<scalarList>("weights"));
        const scalarList abscissae(dict.get<scalarList>("abscissae"));

        if
        (
            weights.size() != abscissae.size()
         || mag(sum(weights) - 1) > 1e-6
         || min(abscissae) < 0
         || max(abscissae) > 1
        )
        {
            FatalIOErrorInFunction(dict)
                << "The nodes of a mixture fraction are as many weights, "
                << "summing to one, as abscissae, within [0, 1]; the "
                << "weights are " << weights << " and the abscissae "
                << abscissae << "." << nl
                << exit(FatalIOError);
        }

        for (label k = 1; k <= maxOrder; k++)
        {
            forAll(weights, nodei)
            {
                m[k] += weights[nodei]/sum(weights)*pow(abscissae[nodei], k);
            }
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "The distribution of the mixture fraction is beta, delta or "
            << "nodes, not " << distribution << "." << nl
            << exit(FatalIOError);
    }

    return m;
}


void Foam::momentGenerationSubModels::mixtureFractionStreams::setMoments
(
    const dictionary& dict,
    const label size
)
{
    // The moments of a composition need one order more of the mixture
    // fraction than their own
    const scalarList m(mixtureFractionMoments(dict, max(xiOrder_) + 1));

    forAll(moments_, mi)
    {
        const label k = xiOrder_[mi];
        const label j = composition_[mi];

        scalar value = m[k];

        if (j == species_.size())
        {
            value = h1_*m[k + 1] + h0_*(m[k] - m[k + 1]);
        }
        else if (j >= 0)
        {
            value = Y1_[j]*m[k + 1] + Y0_[j]*(m[k] - m[k + 1]);
        }

        moments_[mi] = scalarField(size, value);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::mixtureFractionStreams
::mixtureFractionStreams
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const labelListList& momentOrders,
    const label nNodes
)
:
    momentGenerationModel(mesh, dict, phaseName, momentOrders, nNodes),
    species_(readSpecies(dict)),
    Y0_(readStream(dict.subDict("streams").subDict("xi0"))),
    Y1_(readStream(dict.subDict("streams").subDict("xi1"))),
    xiOrder_(nMoments_, 0),
    composition_(nMoments_, -1),
    carriesEnthalpy_(nDimensions_ == species_.size() + 2),
    thermoPtr_(),
    h0_(0),
    h1_(0)
{
    if (species_.size() != nDimensions_ - 1 && !carriesEnthalpy_)
    {
        FatalIOErrorInFunction(dict)
            << "The distribution of " << phaseName << " has "
            << nDimensions_ - 1 << " coordinates conditioned on the "
            << "mixture fraction, and the case has " << species_.size()
            << " species: " << species_ << nl
            << "    The coordinates are the species, or the species and "
            << "the enthalpy after them." << nl
            << exit(FatalIOError);
    }

    if (carriesEnthalpy_)
    {
        readStreamEnthalpies(dict.subDict("streams"));
    }

    // The mixture fraction is the first dimension. A moment is of the
    // mixture fraction alone, or of first order in one composition
    forAll(momentOrders_, mi)
    {
        const labelList& order = momentOrders_[mi];

        xiOrder_[mi] = order[0];

        for (label dimi = 1; dimi < nDimensions_; dimi++)
        {
            if (order[dimi] == 1 && composition_[mi] < 0)
            {
                composition_[mi] = dimi - 1;
            }
            else if (order[dimi] != 0)
            {
                FatalIOErrorInFunction(dict)
                    << "The moment of orders " << order << " is neither of "
                    << "the mixture fraction alone nor of first order in "
                    << "one composition, and compositions linear in the "
                    << "mixture fraction give no rule for it." << nl
                    << exit(FatalIOError);
            }
        }

        // A moment of the enthalpy is an energy per unit mass
        if (composition_[mi] == species_.size())
        {
            momentDimensions_.set
            (
                order,
                new dimensionSet(dimEnergy/dimMass)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationSubModels::mixtureFractionStreams
::~mixtureFractionStreams()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationSubModels::mixtureFractionStreams::updateMoments
(
    const dictionary& dict,
    const label patchi
)
{
    setMoments(dict, reset(patchi));
}


void Foam::momentGenerationSubModels::mixtureFractionStreams::updateMoments
(
    const dictionary& dict,
    const labelList& cells
)
{
    setMoments(dict, reset(cells));
}


void Foam::momentGenerationSubModels::mixtureFractionStreams
::writeDerivedFields
(
    const mappedPtrList<volScalarField>& moments
) const
{
    // The thermo, where it was built, holds the fields of the species
    // under their names: the means are written apart from them
    const volScalarField& m0 = moments(labelList(nDimensions_, 0));

    forAll(species_, j)
    {
        labelList order(nDimensions_, 0);
        order[1 + j] = 1;

        // The mean of the species, with the patches of its moment
        volScalarField Y
        (
            IOobject
            (
                species_[j],
                "0",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            moments(order)
        );

        Y.primitiveFieldRef() /= m0.primitiveField();

        forAll(Y.boundaryField(), patchi)
        {
            Y.boundaryFieldRef()[patchi] ==
                Y.boundaryField()[patchi]/m0.boundaryField()[patchi];
        }

        Info<< "Writing the mean mass fraction of " << species_[j] << endl;

        Y.write();
    }
}


// ************************************************************************* //
