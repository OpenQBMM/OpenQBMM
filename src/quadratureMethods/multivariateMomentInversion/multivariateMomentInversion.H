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

Class
    Foam::multivariateMomentInversion

Description
    Abstract class to allow run-time selection of the algorithm to be used
    for the extended quadrature method of moments.

    References
    \verbatim
        "Conditional quadrature method of moments for kinetic equations"
        C Yuan, R O Fox
        Journal of Computational Physics
        Volume 230, Pages 8216-8246, 2011
    \endverbatim

SourceFiles
    multivariateMomentInversion.C
    multivariateMomentInversionI.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateMomentInversion_H
#define multivariateMomentInversion_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "scalarMatrices.H"
#include "dictionary.H"
#include "univariateMomentSet.H"
#include "Vandermonde.H"
#include "nDimensionalMappedList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class multivariateMomentInversion Declaration
\*---------------------------------------------------------------------------*/
class multivariateMomentInversion
{
    // Private member data

        //- Number of moments
        const label nMoments_;

        //- List of the number of nodes
        const labelList nNodes_;

        //- Number of dimensions
        const label nDims_;

        //- Type of support for each direction
        const List<word> support_;

        //- Quadratrure abscissae
        PtrList< nDimensionalMappedList<scalar> > abscissae_;

        //- Quadratrure weights
        PtrList< nDimensionalMappedList<scalar> > weights_;

        //- Stored moments
        nDimensionalMappedList<scalar> moments_;

        //- Stored conditional moments
        PtrList<PtrList<nDimensionalMappedList<scalar> > > conditionalMoments_;

        //- Inverted VanderMonde matricies
        PtrList< nDimensionalMappedList<scalarSquareMatrix> > invVR_;

    // Private member functions

        //- Cycles through all alphas and solves the given conditional moments
        void cycleAlphaCM
        (
            const label dimi,
            const label dimJ,
            label ai,
            labelList& pos
        );

        //- Sets the Vandermonde matricies for the previous dimension
        void setVR
        (
            const label dimi,
            labelList& pos,
            label ai
        );

        //- Once all of the conditional moments are known, said moments are
        //  inverted to find the weights and abscissae
        void cycleAlphaWheeler
        (
            const label dimi,
            label alphai,
            labelList& pos
        );

        //- Reset inverter
        void reset();

        //- Disallow default bitwise copy construct
        //multivariateMomentInversion(const multivariateMomentInversion&);

        //- Disallow default bitwise assignment
        //void operator=(const multivariateMomentInversion&);


public:

    // Constructors

        // Construct from nMoments, momentMap, nNodes, and support type
        multivariateMomentInversion
        (
            const label nMoments,
            const Map<label> map,
            const labelList& nNodes,
            const List<word>& support
        );


    //- Destructor
    virtual ~multivariateMomentInversion();


    // Member Functions

        //- Invert moments to find weights and  abscissae
        void invert(const nDimensionalMappedList<scalar>& moments);

        //- Return number of moments
        inline label nMoments();

        //- Return number of number of nodes in each dimension
        inline const labelList& nNodes();

        //- Return mapped quadrature weights for each dimension
        inline const PtrList<nDimensionalMappedList<scalar> >& weights() const;

        //- Return mapped quadrature abscissae for each dimension
        inline const PtrList<nDimensionalMappedList<scalar> >& abscissae() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multivariateMomentInversionI.H"

#endif


// ************************************************************************* //
