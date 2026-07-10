/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR AUTHOR
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    Description
        Model for particle deposition in porous media, accounting for:
        - Brownian diffusion
        - Interception mechanisms
        - Hydrodynamic interactions via Kuwabara function
        
    Reference
        Based on filtration theory for particle deposition
\*---------------------------------------------------------------------------*/

#include "porousMediaFiltration.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace depositionModels
{
    defineTypeNameAndDebug(porousMediaFiltration, 0);

    addToRunTimeSelectionTable
    (
        depositionModel,
        porousMediaFiltration,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::depositionModels::porousMediaFiltration::
porousMediaFiltration
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    depositionModel(dict,mesh),
    zoneName_(),
    cellZoneIDs_(),
    epsilon_(dict.lookup("porosity")),
//    s_(dict.lookup("specificArea")),
//    dG_(cellZoneIDs_.size(), Zero),
    dG_(dict.lookup("objSize")),
    ge_(epsilon_.size(), Zero),
    U_(mesh.lookupObject<volVectorField>("U")),
    T_(mesh.lookupObject<volScalarField>("T")),
    mu_(mesh.lookupObject<volScalarField>("thermo:mu"))
{

    dict.readEntry("cellZone", zoneName_);
    cellZoneIDs_ = mesh_.cellZones().indices(zoneName_);
    
    if (returnReduceAnd(cellZoneIDs_.empty()) && Pstream::master())
    {
        FatalErrorInFunction
            << "Cannot find deposition cellZone " << zoneName_ << endl
            << "Valid zones : "
            << flatOutput(mesh_.cellZones().names()) << nl
            << "Valid groups: "
            << flatOutput(mesh_.cellZones().groupNames()) << nl
            << exit(FatalError);
    }

    for (const label zonei : cellZoneIDs_)
    {
//        dG_[zonei] = 6.0/s_[zonei];
        ge_[zonei] =  epsilon_[zonei]/
            (
                    2.0 - epsilon_[zonei]
                - (9.0/5.0)*pow(1.0 - epsilon_[zonei], 1.0/3.0) 
                + 0.2*sqr(1.0 - epsilon_[zonei])
            );
    }

    Info << "depositionModel  minNumDensity : "<<minM0_ << endl;

}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalanceSubModels::depositionModels::porousMediaFiltration
::~porousMediaFiltration()
{}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * //


Foam::scalar 
Foam::populationBalanceSubModels::depositionModels::porousMediaFiltration::Kd
(
    const scalar& dp,
    const label celli,
    const label zonei
) const
{


    // Calculate base deposition rate coefficient
    // b = -(3/2)(|vp|/DG)((1-ε)/ε)η(vp,dp)
    scalar Uc = mag(U_[celli]) + SMALL;

    scalar baseRate = 1.5*Uc/dG_[zonei]*((1.0 - epsilon_[zonei])/epsilon_[zonei]);

    const scalar pi = constant::mathematical::pi;
    const scalar kB = constant::physicoChemical::k.value();

//    scalar Ddiff = kB/3.0/pi*T_[celli]/mu_[celli]/dp;  
    scalar Ddiff = kB/3.0/pi*T_[celli]/mu_[celli]/dG_[zonei];  
//    scalar pe = dp*Uc/Ddiff;  // Peclet number
    scalar pe = dG_[zonei]*Uc/Ddiff;  // Peclet number
    scalar etaD = 3.5*ge_[zonei]*pow(pe, -2.0/3.0);

    scalar mG = dp/dG_[zonei];  // Size ratio
    scalar s = (3.0 - 2.0*epsilon_[zonei])/epsilon_[zonei];
    scalar etaI = 1.5*sqr(mG)*pow(ge_[zonei], 3.0)/(1.0 + mG)*s;

    scalar efficiency = etaD + etaI - etaD*etaI;
/*
    Info << "celli " << celli << tab << "dp " << dp << endl;
    Info << "Uc " << Uc << tab << "T_[celli] " << T_[celli] << tab << "mu_[celli] " << mu_[celli] << endl;
    Info << "ge_  " << ge_ << endl;
    Info << "Ddiff " << Ddiff << endl;
    Info << "pe   " << pe  << endl;
    Info << "etaD   " << etaD  << endl;
    Info << "etaI   " << etaI  << endl;
*/

    return baseRate*efficiency;
    
}


void Foam::populationBalanceSubModels::depositionModels::porousMediaFiltration::correct
(
    scalarQuadratureApproximation& quadrature
) const
{

    PtrList<volScalarNode>& nodes = quadrature.nodes();
    const volScalarMomentFieldSet& moments = quadrature.moments();
    scalar globalDt = mesh_.time().deltaT().value();
    label sizeIndex = nodes[0].sizeIndex();

    forAll(mesh_.C(), celli)
    {
        for (const label zonei : cellZoneIDs_)
        {

            const cellZone& cZone = mesh_.cellZones()[zonei];

            if(cZone.whichCell(celli) != -1 && moments(0)[celli] > minM0_)
            {
                forAll(nodes, pNodei)
                {
                    volScalarNode& node = nodes[pNodei];
            
                    scalar bAbscissa =
                        max(node.primaryAbscissae()[sizeIndex][celli], scalar(0));
                        
                    scalar dp = node.d(celli, bAbscissa);
                    scalar blockRate(min(Kd(dp, celli, zonei) , 20.0));

                    node.primaryWeight()[celli] *= exp(- blockRate*globalDt);
                }

                quadrature.updateLocalMoments(celli);
                    
            }
        }
    }
}

// ************************************************************************* //
