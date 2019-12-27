/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 Alberto Passalacqua
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

#ifndef sizeCHyQMOMMomentInversions_H
#define sizeCHyQMOMMomentInversions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sizeCHyQMOMMomentInversion.H"
#include "CHyQMOMMomentInversion.H"
#include "CHyQMOMPlusMomentInversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multivariateMomentInversions
{
    typedef multivariateMomentInversions::sizeCHyQMOMBase<CHyQMOM> sizeCHyQMOM;
    defineTemplateTypeNameWithName(sizeCHyQMOM, "sizeCHyQMOM");

    typedef sizeCHyQMOMBase<CHyQMOMPlus> sizeCHyQMOMPlus;
    defineTemplateTypeNameWithName(sizeCHyQMOMPlus, "sizeCHyQMOMPlus");
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
