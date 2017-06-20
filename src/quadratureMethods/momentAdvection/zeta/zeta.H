/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 Alberto Passalacqua
     \\/     M anipulation  |
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

Class
    Foam::zeta

Description
    Second-order univariate moment advection with zeta kinetic scheme.

SourceFiles
    zeta.C
    zetaI.H

\*---------------------------------------------------------------------------*/

#ifndef zeta_H
#define zeta_H

#include "univariateMomentAdvection.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class zeta Declaration
\*---------------------------------------------------------------------------*/

class zeta
:
    public univariateMomentAdvection
{
    // Private data

        //- Number of zeta_k values
        label nZetas_;

        //- List of fields of zeta_k (n fields for n + 1 moments)
        PtrList<volScalarField> zetas_;

        //- List of interpolated zeta_k values (neighbour)
        PtrList<surfaceScalarField> zetasNei_;

        //- List of interpolated nodes (owner)
        PtrList<surfaceScalarField> zetasOwn_;

        //- List of interpolated moments (neighbour)
        PtrList<surfaceScalarField> momentsNei_;

        //- List of interpolated moments (owner)
        PtrList<surfaceScalarField> momentsOwn_;

        //- Face velocity
        const surfaceScalarField& phi_;


    // Private member functions

        //- Compute n values of zeta_k from n + 1 moments
        void computeZetaFields();

        //- Updates reconstructed moments from the corresponding values of zeta
        void updateMomentFieldsFromZetas
        (
            const surfaceScalarField m0f,
            const PtrList<surfaceScalarField>& zetaf,
            PtrList<surfaceScalarField>& mf
        );

        //- Compute n + 1 moments from n values of zeta_k assuming m0 = 1.
        void zetaToMoments(const scalarList& zetaf, scalarList& mf);


protected:

        //- Protected member functions
        void interpolateZetas();


public:

    //- Runtime type information
    TypeName("zeta");


    // Constructors

        //- Construct from univariateMomentSet
        zeta
        (
            const dictionary& dict,
            const univariateQuadratureApproximation& quadrature,
            const surfaceScalarField& phi,
            const word& support
        );


    //- Destructor
    virtual ~zeta();


    // Public member functions

        //- Return the maximum Courant number ensuring moment realizability
        virtual scalar realizableCo();

        //- Update moment advection
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
