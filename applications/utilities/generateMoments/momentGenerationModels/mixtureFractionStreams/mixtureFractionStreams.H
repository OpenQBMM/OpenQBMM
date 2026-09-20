/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Alberto Passalacqua
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
    Foam::momentGenerationSubModels::mixtureFractionStreams

Description
    Moments of a mixture fraction and of the compositions conditioned on
    it, for two streams that have only mixed.

    The mixture fraction xi is the mass fraction of the stream at xi = 1
    in a mixture with the stream at xi = 0. Where the two streams have
    mixed and nothing else has happened, every composition is linear in
    the mixture fraction between its values in the streams,

        Y_j(xi) = xi Y_j1 + (1 - xi) Y_j0,

    so the moments of a composition follow from those of the mixture
    fraction, m_k = <xi^k>, alone:

        <xi^k Y_j> = Y_j1 m_(k+1) + Y_j0 (m_k - m_(k+1)).

    This is the state a case starts from and the one its inlets carry,
    and it is exact for the moment set of a quadrature of conditional
    means, which holds <xi^k> for k < 2n and <xi^k Y_j> for k < n, n being
    the number of nodes: the moment of the highest order it needs, m_n, is
    one the set carries. The model generates those moments and refuses a
    set holding any other, which it has no rule for.

    Where the distribution carries one coordinate more after the
    compositions, it is the absolute enthalpy conditioned on the mixture
    fraction, which mixing alone leaves linear in it as well, between the
    enthalpies of the streams at their temperatures and compositions:

        <xi^k h> = h_1 m_(k+1) + h_0 (m_k - m_(k+1)).

    Those enthalpies need the thermodynamic data of the species, which
    the model reaches by building the thermo of the case; they are taken
    at the mean pressure of the case. The thermo reads the mass fraction
    of every species where it is built, so where a case does not have
    those fields yet the model first writes them with the composition of
    the stream at xi = 0, and they are replaced by the means at the end.
    The moments of the enthalpy are given the dimensions of an energy per
    unit mass, whatever abscissaDimension says.

    The compositions are the species of the case, in the order of its
    species list, as the chemistry reaction model takes them: species of
    the dictionary when it is given, and otherwise the species entry of
    the foamChemistryFile named in thermophysicalProperties. The streams
    are given as that reaction model reads them, so a case can keep them
    in one file and include it in both dictionaries; their temperature is
    read only where the enthalpy is carried. The mass fractions of a stream
    are given by species name, those not given are zero, and they must sum
    to one.

    The distribution of the mixture fraction is given for the internal
    field, and for every patch that fixes its value and every region that
    departs from it, as one of
    - beta: the beta distribution of a mean and a variance, 0 < variance
      < mean (1 - mean), whose moments are m_k = m_(k-1) (a + k - 1)/(a +
      b + k - 1), with a = mean c, b = (1 - mean) c and c = mean (1 -
      mean)/variance - 1;
    - delta: all of the fluid at one value, 0 or 1 for an inlet of a pure
      stream;
    - nodes: weights and abscissae, the weights summing to one.
    Every value is uniform over the cells or the faces it is given for.

    The thermo of a reacting case reads the mass fraction of every
    species, which must sum to one where it starts, and the mixing model
    rewrites it at every step with the mean over the environments. The
    model therefore writes, beside the moments, a field for every species
    holding its mean, the moment <Y_j> over the zero-order moment, with
    the patches of that moment.

    \verbatim
    mixing
    {
        type                mixtureFractionStreams;
        weightDimension     [0 0 0 0 0 0 0];
        abscissaDimension   [0 0 0 0 0 0 0];

        species             (O2 H2O CH4 CO2 N2);    // or from the thermo

        streams
        {
            xi0 { T 1200; Y { O2 0.233; N2 0.767; } }
            xi1 { T 300;  Y { CH4 1; } }
        }

        internal    { distribution beta; mean 0.1; variance 0.04; }
        airInlet    { distribution delta; value 0; }
        fuelInlet   { distribution delta; value 1; }
    }
    \endverbatim

SourceFiles
    mixtureFractionStreams.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureFractionStreams_H
#define mixtureFractionStreams_H

#include "momentGenerationModel.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace momentGenerationSubModels
{

/*---------------------------------------------------------------------------*\
                  Class mixtureFractionStreams Declaration
\*---------------------------------------------------------------------------*/

class mixtureFractionStreams
:
    public momentGenerationModel
{
    // Private data

        //- The species, which are the compositions in their order
        wordList species_;

        //- Composition of the stream at xi = 0
        scalarField Y0_;

        //- Composition of the stream at xi = 1
        scalarField Y1_;

        //- For every moment, the order of the mixture fraction
        labelList xiOrder_;

        //- For every moment, the composition it is first order in, -1 for
        //  a moment of the mixture fraction alone, and the number of
        //  species for a moment of the enthalpy
        labelList composition_;

        //- Whether the distribution carries the enthalpy
        const bool carriesEnthalpy_;

        //- The thermo of the case, built for the data of its species
        autoPtr<fluidThermo> thermoPtr_;

        //- Absolute enthalpy of the stream at xi = 0
        scalar h0_;

        //- Absolute enthalpy of the stream at xi = 1
        scalar h1_;


    // Private member functions

        //- The species of the dictionary, or of the chemistry file the
        //  thermophysical properties of the case name
        wordList readSpecies(const dictionary& dict) const;

        //- The composition of a stream
        scalarField readStream(const dictionary& dict) const;

        //- Build the thermo of the case, writing the mass fractions it
        //  reads where the case does not have them yet, and take the
        //  absolute enthalpies of the streams from it
        void readStreamEnthalpies(const dictionary& streams);

        //- The moments of the mixture fraction up to an order, for the
        //  distribution of a dictionary
        scalarList mixtureFractionMoments
        (
            const dictionary& dict,
            const label maxOrder
        ) const;

        //- Set every moment to its uniform value over a number of values
        void setMoments(const dictionary& dict, const label size);


public:

    //- Runtime type information
    TypeName("mixtureFractionStreams");


    // Constructors

        //- Construct from components
        mixtureFractionStreams
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName,
            const labelListList& momentOrders,
            const label nNodes
        );


    //- Destructor
    virtual ~mixtureFractionStreams();


    // Member Functions

        //- Update the moments of the internal field or of a patch
        virtual void updateMoments
        (
            const dictionary& dict,
            const label patchi = -1
        );

        //- Update the moments of a set of cells
        virtual void updateMoments
        (
            const dictionary& dict,
            const labelList& cells
        );

        //- Write the mean mass fraction of every species, for the thermo
        virtual void writeDerivedFields
        (
            const mappedPtrList<volScalarField>& moments
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace momentGenerationSubModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
