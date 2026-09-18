/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2014-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2025 Alberto Passalacqua
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

#include "univariateMomentInversion.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::label Foam::univariateMomentInversion::nQuadratureNodes
(
    const dictionary& dict,
    const label nMoments
)
{
    // The moments determine a quadrature of nMoments/2 nodes, and each
    // inversion adds the nodes of its own: the abscissa a Gauss-Radau
    // quadrature fixes, the two a Gauss-Lobatto one fixes, and the ones the
    // generalized quadrature extends the recurrence relationship to, which
    // the case chooses.
    //
    // This is asked before the inversion exists, since it is what sizes the
    // list of nodes the inversion writes into, so it is the one place that
    // reads the number of nodes of a type from its name rather than from
    // the object.
    const word type(dict.lookup("univariateMomentInversion"));

    label nAdditionalNodes = 0;

    if (type == "GaussRadau")
    {
        nAdditionalNodes = 1;
    }
    else if (type == "GaussLobatto")
    {
        nAdditionalNodes = 2;
    }
    else if (type == "GQMOM")
    {
        nAdditionalNodes =
            dict.get<label>("nAdditionalQuadraturePoints");
    }

    if (nAdditionalNodes < 0)
    {
        FatalIOErrorInFunction(dict)
            << "The number of additional quadrature points cannot be "
            << "negative." << nl
            << exit(FatalIOError);
    }

    return nMoments/2 + nAdditionalNodes;
}


Foam::autoPtr<Foam::univariateMomentInversion>
Foam::univariateMomentInversion::New
(
    const dictionary& dict,
    const label nMaxNodes
)
{
    word univariateMomentInversionType
    (
        dict.lookup("univariateMomentInversion")
    );

    Info<< "Selecting univariateMomentInversion: "
        << univariateMomentInversionType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(univariateMomentInversionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown univariateMomentInversion type "
            << univariateMomentInversionType << endl << endl
            << "Valid univariateMomentInversion types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, nMaxNodes);
}


// ************************************************************************* //
