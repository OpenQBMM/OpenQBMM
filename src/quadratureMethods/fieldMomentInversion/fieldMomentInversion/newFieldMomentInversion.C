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

#include "fieldMomentInversion.H"
#include "univariateMomentInversion.H"
#include "CHyQMOMMomentInversion.H"
#include "CHyQMOMPlusMomentInversion.H"
#include "conditionalMeansMomentInversion.H"

// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

namespace Foam
{

//- The moments and the nodes of a velocity closure of nDimensions velocity
//  directions
static void velocityMomentSet
(
    const word& type,
    const label nDimensions,
    labelListList& momentOrders,
    labelListList& nodeIndexes
)
{
    if (type == "CHyQMOM")
    {
        momentOrders =
            multivariateMomentInversions::CHyQMOM::getMomentOrders
            (
                nDimensions
            );

        nodeIndexes =
            multivariateMomentInversions::CHyQMOM::getNodeIndexes
            (
                nDimensions
            );

        return;
    }

    if (type == "CHyQMOMPlus")
    {
        momentOrders =
            multivariateMomentInversions::CHyQMOMPlus::getMomentOrders
            (
                nDimensions
            );

        nodeIndexes =
            multivariateMomentInversions::CHyQMOMPlus::getNodeIndexes
            (
                nDimensions
            );

        return;
    }

    FatalErrorInFunction
        << "The moment set of " << type << " is not built from the counts "
        << "of the case." << nl
        << "    List the moments and the nodes of the case in "
        << "quadratureProperties." << nl
        << exit(FatalError);
}

} // End namespace Foam



// * * * * * * * * * * * * * * * * * Names * * * * * * * * * * * * * * * * * //

Foam::word Foam::fieldMomentInversion::fieldInversionType
(
    const dictionary& dict
)
{
    const word type(dict.lookup("fieldMomentInversion"));

    if (type == "basicVelocityFieldMomentInversion")
    {
        static bool warned = false;

        if (!warned)
        {
            WarningInFunction
                << "basicVelocityFieldMomentInversion is now named "
                << "basicMultivariateFieldMomentInversion, and its "
                << "dictionary basicMultivariateMomentInversion." << nl
                << "    Rename them in quadratureProperties, or run "
                << "bin/qbmmConvertQuadratureProperties on the case." << nl
                << endl;

            warned = true;
        }

        return "basicMultivariateFieldMomentInversion";
    }

    return type;
}


const Foam::dictionary& Foam::fieldMomentInversion::multivariateInversionDict
(
    const dictionary& dict
)
{
    if (!dict.found("basicMultivariateMomentInversion"))
    {
        if (dict.found("basicVelocityMomentInversion"))
        {
            static bool warned = false;

            if (!warned)
            {
                WarningInFunction
                    << "The dictionary basicVelocityMomentInversion is now "
                    << "named basicMultivariateMomentInversion." << nl
                    << "    Rename it in quadratureProperties, or run "
                    << "bin/qbmmConvertQuadratureProperties on the case."
                    << nl << endl;

                warned = true;
            }

            return dict.subDict("basicVelocityMomentInversion");
        }
    }

    return dict.subDict("basicMultivariateMomentInversion");
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

bool Foam::fieldMomentInversion::buildMomentSet
(
    const dictionary& dict,
    const word& fieldInversion,
    labelListList& momentOrders,
    labelListList& nodeIndexes
)
{
    if (fieldInversion == "basicFieldMomentInversion")
    {
        // A univariate quadrature carries the moments of order zero up, and
        // as many nodes as its inversion builds from them
        const dictionary& inversionDict = dict.subDict("basicMomentInversion");

        const label nMoments =
            inversionDict.lookupOrDefault<label>("nMoments", -1);

        if (nMoments < 0)
        {
            return false;
        }

        const label nNodes =
            univariateMomentInversion::nQuadratureNodes
            (
                inversionDict,
                nMoments
            );

        momentOrders.setSize(nMoments);

        forAll(momentOrders, mi)
        {
            momentOrders[mi] = labelList(1, mi);
        }

        nodeIndexes.setSize(nNodes);

        forAll(nodeIndexes, nodei)
        {
            nodeIndexes[nodei] = labelList(1, nodei);
        }

        return true;
    }

    if (fieldInversion != "basicMultivariateFieldMomentInversion")
    {
        return false;
    }

    const dictionary& inversionDict = multivariateInversionDict(dict);

    const word type(inversionDict.lookup("type"));

    // Conditional means: the moments of the primary coordinate and the
    // first-order moments of every coordinate conditioned on it, the
    // scalars of the model and the velocities of the case
    if (type == "conditionalMeans")
    {
        const label nMoments =
            inversionDict.lookupOrDefault<label>("nMoments", -1);

        const label nScalars =
            inversionDict.lookupOrDefault<label>("nScalars", -1);

        if (nMoments < 0 || nScalars < 0)
        {
            return false;
        }

        multivariateMomentInversions::conditionalMeans::momentOrdersAndNodes
        (
            nMoments,
            nScalars,
            inversionDict.lookupOrDefault<label>("nVelocityDimensions", 0),
            momentOrders,
            nodeIndexes
        );

        return true;
    }

    const label nVelocityDimensions =
        inversionDict.lookupOrDefault<label>("nVelocityDimensions", -1);

    if (nVelocityDimensions < 0)
    {
        return false;
    }

    // A monokinetic closure carries the moments of the size alone, and the
    // first-order velocity moments times the powers of the size below the
    // number of its nodes: one velocity at every size node
    if (type == "monoKinetic")
    {
        const label nSizeMoments =
            inversionDict.lookupOrDefault<label>("nSizeMoments", -1);

        if (nSizeMoments < 0)
        {
            return false;
        }

        multivariateMomentInversions::conditionalMeans::momentOrdersAndNodes
        (
            nSizeMoments,
            0,
            nVelocityDimensions,
            momentOrders,
            nodeIndexes
        );

        return true;
    }

    // The closure whose moments the case carries. It is the one that
    // inverts them unless the case asks for another, as a case that
    // transports the moments of CHyQMOMPlus and inverts them with CHyQMOM
    // does.
    const word transported
    (
        inversionDict.lookupOrDefault<word>
        (
            "transportedMoments",
            type == "sizeCHyQMOM"
          ? inversionDict.lookupOrDefault<word>("sizeVelocityClosure", "CHyQMOM")
          : type
        )
    );

    labelListList velocityOrders;
    labelListList velocityNodes;

    velocityMomentSet
    (
        transported,
        nVelocityDimensions,
        velocityOrders,
        velocityNodes
    );

    if (type != "sizeCHyQMOM")
    {
        momentOrders = velocityOrders;
        nodeIndexes = velocityNodes;

        return true;
    }

    // A size-velocity quadrature carries the moments of the size alone, and
    // the velocity moments of the closure conditioned on each of the size
    // nodes the size moments build
    const label nSizeMoments =
        inversionDict.lookupOrDefault<label>("nSizeMoments", -1);

    if (nSizeMoments < 0)
    {
        return false;
    }
    const label nSizeNodes = nSizeMoments/2;

    momentOrders.clear();
    nodeIndexes.clear();

    for (label k = 0; k < nSizeMoments; k++)
    {
        labelList order(1 + nVelocityDimensions, Zero);
        order[0] = k;
        momentOrders.append(order);
    }

    for (label k = 0; k < nSizeNodes; k++)
    {
        forAll(velocityOrders, vi)
        {
            bool zero = true;

            forAll(velocityOrders[vi], dimi)
            {
                zero = zero && velocityOrders[vi][dimi] == 0;
            }

            if (zero)
            {
                continue;
            }

            labelList order(1 + nVelocityDimensions, Zero);
            order[0] = k;

            forAll(velocityOrders[vi], dimi)
            {
                order[dimi + 1] = velocityOrders[vi][dimi];
            }

            momentOrders.append(order);
        }
    }

    for (label k = 0; k < nSizeNodes; k++)
    {
        forAll(velocityNodes, vi)
        {
            labelList index(1 + nVelocityDimensions, Zero);
            index[0] = k;

            forAll(velocityNodes[vi], dimi)
            {
                index[dimi + 1] = velocityNodes[vi][dimi];
            }

            nodeIndexes.append(index);
        }
    }

    return true;
}


Foam::List<Foam::supportType> Foam::fieldMomentInversion::momentSupports
(
    const dictionary& dict,
    const label nDimensions,
    const List<supportType>& modelSupports
)
{
    const word fieldInversion(fieldInversionType(dict));

    if (fieldInversion != "basicMultivariateFieldMomentInversion")
    {
        return modelSupports;
    }

    // The velocity models used to give every quadrature three dimensions
    // on the reals, whatever the number of its dimensions and whether one
    // of them is a size
    List<supportType> supports(nDimensions, supportType::R);

    const word type(multivariateInversionDict(dict).lookup("type"));

    // Conditional means are taken on the supports of the model, which is
    // what the closure of a mixture fraction and its compositions needs
    if (type == "conditionalMeans")
    {
        return modelSupports;
    }

    if (type == "sizeCHyQMOM" || type == "monoKinetic")
    {
        supports[0] = supportType::RPlus;
    }

    return supports;
}


void Foam::fieldMomentInversion::momentSet
(
    const dictionary& dict,
    labelListList& momentOrders,
    labelListList& nodeIndexes
)
{
    const word fieldInversion(fieldInversionType(dict));

    // A case written before the inversions owned their moment sets lists
    // them itself. It is taken as it is, so that nothing it carries is
    // dropped, and it is told what its inversion would build for it.
    if (dict.found("moments"))
    {
        momentOrders = labelListList(dict.lookup("moments"));
        nodeIndexes = labelListList(dict.lookup("nodes"));

        labelListList builtOrders;
        labelListList builtNodes;

        if (buildMomentSet(dict, fieldInversion, builtOrders, builtNodes))
        {
            if (builtOrders == momentOrders && builtNodes == nodeIndexes)
            {
                WarningInFunction
                    << "The moments and the nodes of this quadrature are "
                    << "the ones its inversion builds, so the lists of them "
                    << "are no longer needed." << nl
                    << "    Remove them from quadratureProperties." << nl
                    << endl;
            }
            else
            {
                WarningInFunction
                    << "This quadrature lists " << momentOrders.size()
                    << " moments and " << nodeIndexes.size()
                    << " nodes, where its inversion builds "
                    << builtOrders.size() << " and " << builtNodes.size()
                    << "." << nl
                    << "    The lists are used. A case that carries the "
                    << "moments of another closure says so with "
                    << "transportedMoments." << nl
                    << endl;
            }
        }

        return;
    }

    if (!buildMomentSet(dict, fieldInversion, momentOrders, nodeIndexes))
    {
        FatalIOErrorInFunction(dict)
            << "The moment set of " << fieldInversion << " is not built "
            << "from the counts of the case." << nl
            << "    List the moments and the nodes of the case in "
            << "quadratureProperties." << nl
            << exit(FatalIOError);
    }
}


Foam::autoPtr<Foam::fieldMomentInversion> Foam::fieldMomentInversion::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelListList& momentOrders,
    const labelListList& nodeIndexes,
    const labelList& velocityIndexes
)
{
    const word fieldMomentInversionType(fieldInversionType(dict));

    Info<< "Selecting fieldMomentInversion: "
        << fieldMomentInversionType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(fieldMomentInversionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fieldMomentInversion type "
            << fieldMomentInversionType << endl << endl
            << "Valid fieldMomentInversion types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()
    (
        dict,
        mesh,
        momentOrders,
        nodeIndexes,
        velocityIndexes
    );
}


// ************************************************************************* //
