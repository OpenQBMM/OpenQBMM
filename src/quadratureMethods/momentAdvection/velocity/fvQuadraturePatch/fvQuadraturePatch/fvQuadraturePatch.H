/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 Alberto Passalacqua
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
    Foam::fvQuadraturePatch

Description
    Abstract class for the advection of a velocity moment set.

SourceFiles
    fvQuadraturePatch.C
    newFvQuadraturePatch.C

\*---------------------------------------------------------------------------*/

#ifndef fvQuadraturePatch_H
#define fvQuadraturePatch_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "quadratureNodes.H"
#include "momentFieldSets.H"
#include "quadratureApproximations.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class fvQuadraturePatch Declaration
\*---------------------------------------------------------------------------*/

class fvQuadraturePatch
{
    // Private member functions

        //- Disallow default bitwise copy construct
        fvQuadraturePatch(const fvQuadraturePatch&);

        //- Disallow default bitwise assignment
        void operator=(const fvQuadraturePatch&);


protected:

    // Protected data

        //- Name of the distribution associated to the quadrature approximation.
        const fvPatch& patch_;

        //- Patch index
        const label patchi_;

        //- Reference to the moments to advect
        const velocityQuadratureApproximation& quadrature_;

        //- Owner nodes
        PtrList<surfaceVelocityNode>& nodesOwn_;

        //- Neighbour nodes
        PtrList<surfaceVelocityNode>& nodesNei_;


    // Protected member functions

        //- Copy the weight of the cell next to each face of the patch onto
        //  both sides of that face. Every condition that reflects or
        //  filters the abscissae leaves the weights alone to begin with,
        //  and the two that rebalance the flux do so afterwards.
        //
        //  The scalar abscissae are deliberately not touched here or in
        //  any condition. They are interpolated onto the faces by the
        //  advection scheme, and reflecting a node off a wall does not
        //  change the size of what it carries.
        void copyWeightsFromCells();


public:

    //- Runtime type information
    TypeName("fvQuadraturePatch");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvQuadraturePatch,
            dictionary,
            (
                const fvPatch& patch,
                const dictionary& dict,
                const velocityQuadratureApproximation& quadrature,
                PtrList<surfaceVelocityNode>& nodesOwn,
                PtrList<surfaceVelocityNode>& nodesNei
            ),
            (patch, dict, quadrature, nodesOwn, nodesNei)
        );


    // Constructors

        //- Construct from velocityMomentSet
        fvQuadraturePatch
        (
            const fvPatch& patch,
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            PtrList<surfaceVelocityNode>& nodesOwn,
            PtrList<surfaceVelocityNode>& nodesNei
        );


    //- Destructor
    virtual ~fvQuadraturePatch();


    // Selectors

        static autoPtr<fvQuadraturePatch> New
        (
            const fvPatch& patch,
            const dictionary& dict,
            const velocityQuadratureApproximation& quadrature,
            PtrList<surfaceVelocityNode>& nodesOwn,
            PtrList<surfaceVelocityNode>& nodesNei
        );


    // Member Functions

        //- Update boundary flux
        virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
