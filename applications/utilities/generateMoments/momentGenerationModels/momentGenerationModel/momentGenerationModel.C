/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 Alberto Passalacqua
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
        moments_[mi].value() = 0.0;

        for (label nodei = 0; nodei < nNodes_; nodei++)
        {
            scalar absCmpt = 1.0;
            forAll(abscissae_[0], cmpti)
            {
                absCmpt *=
                    pow
                    (
                        abscissae_[nodei][cmpti],
                        momentOrder[cmpti]
                    );
            }
            moments_[mi].value() += weights_[nodei]*absCmpt;
        }
    }
}

Foam::labelListList Foam::momentGenerationModel::createNodeIndexes() const
{
    labelListList nodeIndexes(0);
    if (nDims_ == 1)
    {
        nodeIndexes = labelListList(nNodes_, labelList(1));
        forAll(nodeIndexes, nodei)
        {
            nodeIndexes[nodei][0] = nodei;
        }
    }
    else if (nDims_ == 2 && nNodes_ == 9)
    {
        nodeIndexes =
        {
            {1, 1},
            {1, 2},
            {1, 3},
            {2, 1},
            {2, 2},
            {2, 3},
            {3, 1},
            {3, 2},
            {3, 3}
        };
    }
    else if (nDims_ == 3 && nNodes_ == 27)
    {
        nodeIndexes =
        {
            {1, 1, 1},
            {1, 1, 2},
            {1, 1, 3},
            {1, 2, 1},
            {1, 2, 2},
            {1, 2, 3},
            {1, 3, 1},
            {1, 3, 2},
            {1, 3, 3},
            {2, 1, 1},
            {2, 1, 2},
            {2, 1, 3},
            {2, 2, 1},
            {2, 2, 2},
            {2, 2, 3},
            {2, 3, 1},
            {2, 3, 2},
            {2, 3, 3},
            {3, 1, 1},
            {3, 1, 2},
            {3, 1, 3},
            {3, 2, 1},
            {3, 2, 2},
            {3, 2, 3},
            {3, 3, 1},
            {3, 3, 2},
            {3, 3, 3}
        };
    }
    else
    {
        FatalErrorInFunction
            << nDims_ << " dimensions and " << nNodes_
            << " nodes does not have a default set of indicies."
            << abort(FatalError);
    }
    return nodeIndexes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::momentGenerationModel::momentGenerationModel
(
    const dictionary& dict,
    const labelListList& momentOrders,
    const label nNodes
)
:
    dict_(dict),
    nDims_(momentOrders[0].size()),
    nNodes_(nNodes),
    nMoments_(momentOrders.size()),
    nodeIndexes_
    (
        dict.found("nodeIndexes")
      ? dict.lookup("nodeIndexes")
      : createNodeIndexes()
    ),
    momentOrders_(momentOrders),
    weights_(nNodes_),
    abscissae_(nNodes_, scalarList(nDims_, 0.0)),
    moments_(nMoments_, momentOrders_)
{
    dimensionSet wDims(dict.lookup("weightDimension"));
    forAll(moments_, mi)
    {
        const labelList& momentOrder = momentOrders_[mi];
        dimensionSet absCmptDims(dimless);
        forAll(momentOrders_[mi], cmpti)
        {
            word absName ="abscissaeDim" + Foam::name(cmpti) + "Dimensions";
            dimensionSet absDim
            (
                (
                    dict.found(absName)
                  ? dict.lookup(absName)
                  : dict.lookup("abscissaDimension")
                )
            );
            absCmptDims *= pow(absDim, momentOrder[cmpti]);
        }
        moments_[mi].dimensions().reset(wDims*absCmptDims);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::momentGenerationModel::~momentGenerationModel()
{}


// ************************************************************************* //
