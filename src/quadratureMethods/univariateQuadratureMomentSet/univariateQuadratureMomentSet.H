/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 Alberto Passalacqua
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
    Foam::univariateQuadratureMomentSet

Description
    Stores a moment set and the associated quadrature.

SourceFiles
    univariateQuadratureMomentSet.C
    univariateQuadratureMomentSet.H

\*---------------------------------------------------------------------------*/

#ifndef univariateQuadratureMomentSet_H
#define univariateQuadratureMomentSet_H

#include "univariateMomentSet.H"
#include "univariateMomentInversion.H"
#include "dictionary.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class univariateMomentInversion;

/*---------------------------------------------------------------------------*\
                    Class univariateQuadratureMomentSet Declaration
\*---------------------------------------------------------------------------*/

class univariateQuadratureMomentSet
:
    public univariateMomentSet
{
    // Private data

        //- Quadrature abscissae
        scalarList abscissae_;

        //- Quadrature weights
        scalarList weights_;

        //- Bool to indicate the moment set was inverted
        bool inverted_;

        //- Minimum known abscissa
        scalar minKnownAbscissa_;

        //- Maximum known abscissa
        scalar maxKnownAbscissa_;

        //- Moment inversion method
        autoPtr<univariateMomentInversion> momentInverter_;


public:

    // Constructors

        //- Construct from number of moments and initial value
        univariateQuadratureMomentSet
        (
            const dictionary& dict,
            const label nMoments,
            const word& support,
            const label nFixedQuadraturePoints = 0,
            const scalar minKnownAbscissa = 0,
            const scalar maxKnownAbscissa = 0
        );

        //- Construct from scalarList
        univariateQuadratureMomentSet
        (
            const dictionary& dict,
            const scalarList& m,
            const word& support,
            const label nFixedQuadraturePoints = 0,
            const scalar minKnownAbscissa = 0,
            const scalar maxKnownAbscissa = 0
        );

    //- Destructor
    virtual ~univariateQuadratureMomentSet();


    // Member Functions

        //- Return quadrature abscissae
        inline const scalarList& abscissae() const;

        //- Compute the quadrature associated to the moment set
        void invert();

        //- Update moments from the corresponding quadrature
        virtual void update();

        //- Return quadrature abscissae
        inline const scalarList& weights() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "univariateQuadratureMomentSetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
