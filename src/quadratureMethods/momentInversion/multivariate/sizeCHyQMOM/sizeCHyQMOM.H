/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 Alberto Passalacqua
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
    Foam::sizeCHyQMOM

Description
    Implementation of the hyperbolic conditional quadrature method of moment
    for velocity with conditioning on size.

SourceFiles
    sizeCHyQMOM.C
    sizeCHyQMOMI.H

\*---------------------------------------------------------------------------*/

#ifndef sizeCHyQMOM_H
#define sizeCHyQMOM_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "vector2D.H"
#include "vector.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "mappedList.H"
#include "Vandermonde.H"
#include "hyperbolicConditionalMomentInversion.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class sizeCHyQMOM Declaration
\*---------------------------------------------------------------------------*/
class sizeCHyQMOM
{
    // Private member data

        //- Number of dimensions in distribution
        const label nDistributionDims_;

        //- Number of dimensions in velocity space
        const label nGeometricDimensions_;

        //- Number of moments
        const label nMoments_;

        //- Number of size moments
        const label nSizeMoments_;

        //- Number of velocity moments
        const label nVelocityMoments_;

        //- Moment orders of full distribution
        const labelListList& momentOrders_;

        //- Moment orders of conditional velocity distribution
        const labelListList velocityMomentOrders_;

        //- Total number of nodes
        const label nNodes_;

        //-Number of nodes in size direction
        const label nSizeNodes_;

        //- Node indexes of full distribution
        const labelListList& nodeIndexes_;

        //- Node indexes of conditional velocity distribution
        const labelListList velocityNodeIndexes_;

        //- Support of the NDF
        const wordList supports_;

        //- Quadratrure weights
        mappedList<scalar> weights_;

        //- Quadratrure abscissae
        mappedList<scalarList> abscissae_;

        //- Univariate moment inversion method
        autoPtr<univariateMomentInversion> sizeInverter_;

        //- Velocity moment inversion method
        autoPtr<hyperbolicConditionalMomentInversion> velocityInverter_;


    // Private member functions

        //- Calculate the number of pure size moments
        label calcNSizeMoments(const labelListList& momentOrders);

        //- Reset inverter
        void reset();

        //- Disallow default bitwise copy construct
        //hyperbolicConditionalMomentInversion(const hyperbolicConditionalMomentInversion&);

        //- Disallow default bitwise assignment
        //void operator=(const hyperbolicConditionalMomentInversion&);


public:

    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        sizeCHyQMOM
        (
            const dictionary& dict,
            const labelListList& momentOrders,
            const labelListList& nodeIndexes
        );


    //- Destructor
    virtual ~sizeCHyQMOM();


    // Member Functions

        //- Invert moments to find weights and  abscissae
        void invert(const multivariateMomentSet& moments);

        //- Return number of moments
        inline label nMoments();

        //- Return number of number of nodes in each dimension
        inline label nNodes();

        //- Return mapped quadrature weights for each dimension
        inline const mappedList<scalar>& weights() const;

        //- Return mapped quadrature abscissae for each dimension
        inline const mappedList<scalarList>& abscissae() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sizeCHyQMOMI.H"

#endif


// ************************************************************************* //
