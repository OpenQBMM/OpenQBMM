/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 Alberto Passalacqua
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::PtrList<Foam::extendedVolScalarNode>&
Foam::univariateQuadratureApproximation::nodes() const
{
    return nodes_();
}

const Foam::momentFieldSet<Foam::volUnivariateMoment, Foam::extendedVolScalarNode>& 
Foam::univariateQuadratureApproximation::moments() const
{
    return moments_;
}

Foam::momentFieldSet<Foam::volUnivariateMoment, Foam::extendedVolScalarNode>& 
Foam::univariateQuadratureApproximation::moments()
{
    return moments_;
}

const Foam::momentFieldSet
<
    Foam::surfaceUnivariateMoment, Foam::extendedSurfaceScalarNode
>& 
Foam::univariateQuadratureApproximation::momentsNei() const
{
    return momentsNei_;
}

Foam::momentFieldSet
<
    Foam::surfaceUnivariateMoment, Foam::extendedSurfaceScalarNode
>&  
Foam::univariateQuadratureApproximation::momentsNei()
{
    return momentsNei_;
}

const Foam::momentFieldSet
<
    Foam::surfaceUnivariateMoment, Foam::extendedSurfaceScalarNode
>& 
Foam::univariateQuadratureApproximation::momentsOwn() const
{
    return momentsOwn_;
}

Foam::momentFieldSet
<
    Foam::surfaceUnivariateMoment, Foam::extendedSurfaceScalarNode
>&  
Foam::univariateQuadratureApproximation::momentsOwn()
{
    return momentsOwn_;
}

Foam::label Foam::univariateQuadratureApproximation::nDimensions() const
{
    return nDimensions_;
}

Foam::label Foam::univariateQuadratureApproximation::nMoments() const
{
    return nMoments_;
}

// ************************************************************************* //
