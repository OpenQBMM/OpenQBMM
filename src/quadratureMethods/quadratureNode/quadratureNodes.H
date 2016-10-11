/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 Alberto Passalacqua
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

InClass
    Foam::quadratureNodes

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef quadratureNodes_H
#define quadratureNodes_H

#include "volFields.H"
#include "quadratureNode.H"
#include "extendedQuadratureNode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef quadratureNode<volScalarField, volScalarField> basicVolScalarNode;

typedef quadratureNode<surfaceScalarField, surfaceScalarField>
    basicSurfaceScalarNode;

typedef quadratureNode<volScalarField, volVectorField> basicVolVectorNode;

typedef quadratureNode<surfaceScalarField, surfaceVectorField>
    basicSurfaceVectorNode;

typedef extendedQuadratureNode<volScalarField, volScalarField, volScalarField>
    extendedVolScalarNode;

typedef extendedQuadratureNode
<
    surfaceScalarField, surfaceScalarField, surfaceScalarField
>
extendedSurfaceScalarNode;

typedef extendedQuadratureNode<volScalarField, volVectorField, volScalarField>
    extendedVolVectorNode;

typedef extendedQuadratureNode
<
    surfaceScalarField, surfaceVectorField, surfaceScalarField
>
extendedSurfaceVectorNode;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
