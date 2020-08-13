/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "phaseCompressibleMeanVelocityForce.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(phaseCompressibleMeanVelocityForce, 0);

    addToRunTimeSelectionTable
    (
        option,
        phaseCompressibleMeanVelocityForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::phaseCompressibleMeanVelocityForce::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::phaseCompressibleMeanVelocityForce::phaseCompressibleMeanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    alpha_
    (
        mesh.thisDb().lookupObject<volScalarField>
        (
            coeffs_.lookupOrDefault<word>("alphaName", "alpha")
        )
    ),
    rho_
    (
        mesh.thisDb().lookupObject<volScalarField>
        (
            coeffs_.lookupOrDefault<word>("rhoName", "rho")
        )
    ),
    Ubar_(coeffs_.lookup("Ubar")),
    magUbar_(max(mag(Ubar_), SMALL) ),
    flowDir_(Ubar_/magUbar_),
    relaxation_(coeffs_.lookupOrDefault<scalar>("relaxation", 1.0)),
    gradP0_(0.0),
    dGradP_(0.0),
    rAPtr_(nullptr)
{
    coeffs_.lookup("fieldNames") >> fieldNames_;

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP0_;
    }

    Info<< "    Initial pressure gradient = " << gradP0_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::phaseCompressibleMeanVelocityForce::correct(volVectorField& U)
{
    const scalarField& rAU = rAPtr_().internalField();

    // Integrate flow variables over cell set
    scalar magAlphaRhoUAve = 0.0;
    scalar alphaRhoAve = 0.0;
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        scalar volCell = cv[cellI];

        magAlphaRhoUAve +=
            rho_[cellI]*alpha_[cellI]*(flowDir_ & U[cellI])*volCell;

        alphaRhoAve += rho_[cellI]*alpha_[cellI]*volCell;
        rAUave += rAU[cellI]*volCell;
    }

    // Collect across all processors
    reduce(magAlphaRhoUAve, sumOp<scalar>());
    reduce(alphaRhoAve, sumOp<scalar>());
    reduce(rAUave, sumOp<scalar>());

    // Volume averages
    magAlphaRhoUAve /= V_;
    alphaRhoAve /= V_;
    rAUave /= V_;

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    dGradP_ = relaxation_*(alphaRhoAve*magUbar_ - magAlphaRhoUAve)/rAUave;

    // Apply correction to velocity field
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        U[cellI] += flowDir_*rAU[cellI]*dGradP_;
    }

    scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected Mean Velocity Magnitude = "
        << magAlphaRhoUAve/alphaRhoAve
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);
}


void Foam::fv::phaseCompressibleMeanVelocityForce::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, vector::zero)
    );

    scalar gradP = gradP0_ + dGradP_;

    UIndirectList<vector>(Su, cells_) = flowDir_*gradP;

    eqn += Su;
}


void Foam::fv::phaseCompressibleMeanVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
    if (!rAPtr_.valid())
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":invA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                1.0/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1.0/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0.0;
}


// ************************************************************************* //
