/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 Alberto Passalacqua
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2021 Alberto Passalacqua
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

#include "momentGenerationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(momentGenerationModel, 0);
    defineRunTimeSelectionTable(momentGenerationModel, dictionary);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::momentGenerationModel::updateMoments()
{
    forAll(moments_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        moments_[mi] = 0.0;

        forAll(abscissae_, nodei)
        {
            scalarField mCmpt = weights_[nodei];

            forAll(abscissae_[0], cmpti)
            {
                mCmpt *=
                    pow
                    (
                        abscissae_[nodei][cmpti],
                        momentOrder[cmpti]
                    );
            }

            moments_[mi] += mCmpt;
        }
    }
}

Foam::label Foam::momentGenerationModel::reset(const label patchi)
{
    label size =
    (
        patchi == -1 ? mesh_.nCells() : mesh_.boundaryMesh()[patchi].size()
    );

    forAll(abscissae_, nodei)
    {
        forAll(abscissae_[nodei], cmpti)
        {
            abscissae_[nodei][cmpti] = scalarField(size, Zero);
        }

        weights_[nodei] = scalarField(size, Zero);
    }

    forAll(moments_, mi)
    {
        moments_[mi] = scalarField(size, Zero);
    }

    return size;
}


Foam::label Foam::momentGenerationModel::reset(const labelList& cells)
{
    label size = cells.size();

    forAll(abscissae_, nodei)
    {
        forAll(abscissae_[nodei], cmpti)
        {
            abscissae_[nodei][cmpti] = scalarField(size, Zero);
        }

        weights_[nodei] = scalarField(size, Zero);
    }

    forAll(moments_, mi)
    {
        moments_[mi] = scalarField(size, Zero);
    }

    return size;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationModel::momentGenerationModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const labelListList& momentOrders,
    const label nNodes
)
:
    mesh_(mesh),
    dict_(dict),
    nDimensions_(momentOrders[0].size()),
    nNodes_(nNodes),
    nMoments_(momentOrders.size()),
    momentOrders_(momentOrders),
    weights_(nNodes_),
    abscissae_(nNodes_, List<scalarField>(nDimensions_)),
    momentDimensions_(nMoments_, momentOrders_),
    moments_(nMoments_, momentOrders_)
{
    dimensionSet weightDimension_(dict.lookup("weightDimension"));
    
    forAll(moments_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        dimensionSet abscissaCmptDimensions(dimless);

        forAll(momentOrders_[mi], cmpti)
        {
            word abscissaName ="abscissaDimension" + Foam::name(cmpti);
        
            dimensionSet abscissaDimension_
            (
                (
                    dict.found(abscissaName)
                  ? dict.lookup(abscissaName)
                  : dict.lookup("abscissaDimension")
                )
            );
            
            abscissaCmptDimensions *= 
                pow(abscissaDimension_, momentOrder[cmpti]);
        }
        momentDimensions_.set
        (
            momentOrder, 
            new dimensionSet(weightDimension_*abscissaCmptDimensions)
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationModel::~momentGenerationModel()
{}


// ************************************************************************* //
