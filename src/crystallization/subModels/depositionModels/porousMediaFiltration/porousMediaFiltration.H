/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR AUTHOR
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
\*---------------------------------------------------------------------------*/


#ifndef porousMediaFiltration_H
#define porousMediaFiltration_H


#include "depositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalanceSubModels
{
namespace depositionModels
{

/*---------------------------------------------------------------------------*\
                    Class porousMediaFiltration Declaration
\*---------------------------------------------------------------------------*/

class porousMediaFiltration
:
    public depositionModel
{

private:


        //- Name(s) of cell-zone
        wordRe zoneName_;

        //- Cell zone IDs
        labelList cellZoneIDs_;

        // Private data
        const scalarList epsilon_;      // porosity

//        const scalarList s_;            //specific surface area

        const scalarList dG_;          // average diameter

        scalarList ge_;          // Kuwabara function

        const volVectorField&  U_;

        //- Temperature field
        const volScalarField&  T_;

        const volScalarField&  mu_;


        scalar Kd
        (
            const scalar& dp,
            const label celli,
            const label zonei
        ) const;

public:

        //- Runtime type information
        TypeName("porousMediaFiltration");

    // Constructors

        //- Construct from components

        porousMediaFiltration
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        // Member Functions
        //- Destructor
        virtual ~porousMediaFiltration();


    // Member Functions

        //- Return Breakup source term
        virtual void correct
        (
            scalarQuadratureApproximation& quadrature
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupKernels
} // End namespace populationBalanceSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //