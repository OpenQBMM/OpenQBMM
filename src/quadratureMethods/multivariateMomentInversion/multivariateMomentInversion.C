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

\*---------------------------------------------------------------------------*/

#include "multivariateMomentInversion.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multivariateMomentInversion::multivariateMomentInversion
(
    const label nMoments,
    const Map<label> map,
    const labelList& nNodes,
    const List<word>& support
)
:
    nMoments_(nMoments),
    nNodes_(nNodes),
    nDims_(nNodes_.size()),
    support_(support),
    abscissae_(nDims_),
    weights_(nDims_),
    moments_(nMoments_, nDims_, map),
    conditionalMoments_(nDims_),
    invVR_(nDims_ - 1)
{
    labelList nNodesCM = nNodes_;

    forAll(nNodes_, dimi)
    {
        weights_.set
        (
            dimi,
            new nDimensionalMappedList<scalar>
            (
                dimi + 1,
                nNodes_
            )
        );

        abscissae_.set
        (
            dimi,
            new nDimensionalMappedList<scalar>
            (
                dimi + 1,
                nNodes_
            )
        );

        forAll(abscissae_[dimi], ai)
        {
            weights_[dimi].set
            (
                ai,
                new scalar(0.0)
            );

            abscissae_[dimi].set
            (
                ai,
                new scalar(0.0)
            );
        }

    }

    forAll(conditionalMoments_, dimi)
    {
        conditionalMoments_.set
        (
            dimi,
            new PtrList<nDimensionalMappedList<scalar> >
            (
                dimi + 1
            )
        );

        nNodesCM = nNodes_;
        nNodesCM[dimi] = 2*nNodes_[dimi];

        forAll(conditionalMoments_[dimi], dimj)
        {
            conditionalMoments_[dimi].set
            (
                dimj,
                new nDimensionalMappedList<scalar>
                (
                    dimi + 1,
                    nNodesCM
                )
            );

            forAll(conditionalMoments_[dimi][dimj], ai)
            {
                conditionalMoments_[dimi][dimj].set
                (
                    ai,
                    new scalar(0.0)
                );
            }
        }
    }
    forAll(invVR_, dimi)
    {
        invVR_.set
        (
            dimi,
            new nDimensionalMappedList<scalarSquareMatrix>
            (
                dimi,
                nNodes_
            )
        );

        forAll(invVR_[dimi], ai)
        {
            invVR_[dimi].set
            (
                ai,
                new scalarSquareMatrix
                (
                    nNodes_[dimi],
                    scalar(0.0)
                )
            );
        }
    }

    forAll(moments_, mi)
    {
        moments_.set
        (
            mi,
            new scalar(0.0)
        );
    }

    // Set all lists to zero
    reset();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multivariateMomentInversion::~multivariateMomentInversion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multivariateMomentInversion
::invert(const nDimensionalMappedList<scalar>& moments)
{
    reset();

    forAll(moments, mi)
    {
        moments_[mi] = moments[mi];
    }

    // Invert primary direction first
    labelList pos(nDims_, 0);
    univariateMomentSet momentsToInvert
    (
        2*nNodes_[0],
        0.0,
        "Gauss",
        support_[0]
    );

    for (label mi = 0; mi < 2*nNodes_[0]; mi++)
    {
        pos[0] = mi;
        momentsToInvert[mi] = moments_(pos);
    }

    momentsToInvert.invert();

    for (label nodei = 0; nodei < nNodes_[0]; nodei++)
    {
        // First component of the weights and abscissae only have one direction
        // Copy from univariateMomentSet to stored copy
        weights_[0](nodei) = momentsToInvert.weights()[nodei];
        abscissae_[0](nodei) = momentsToInvert.abscissae()[nodei];
    }

    // Solve remaining directions
    for (label dimi = 1; dimi < nDims_; dimi++)
    {
        for (label dimj = 0; dimj < dimi; dimj++)
        {
            if (dimj == 0)
            {
                labelList pos(dimi, 0);
                setVR(dimi - 1, pos, 0);
            }
            labelList posC(dimi + 1, 0);
            cycleAlphaCM(dimi, dimj, 0, posC);
        }

        labelList posW(dimi + 1, 0);
        cycleAlphaWheeler(dimi, 0, posW);
    }
}

void Foam::multivariateMomentInversion::cycleAlphaCM
(
    const label dimi,
    const label dimj,
    label ai,
    labelList& pos
)
{
    if (dimj == ai)
    {
        cycleAlphaCM(dimi,dimj,ai+1,pos);
        return;
    }
    else if (ai < dimi)
    {
        for (label i = 0; i < nNodes_[ai]; i++)
        {
            pos[ai] = i;
            cycleAlphaCM(dimi,dimj,ai+1,pos);
        }
        return;
    }
    else if (ai == dimi)
    {
        for (label i = 0; i < 2*nNodes_[dimi]; i++)
        {
            pos[dimi] = i;
            cycleAlphaCM(dimi,dimj,ai+1,pos);
        }
        return;
    }
    else
    {
        scalarRectangularMatrix Yold(nNodes_[dimj], 1, 0.0);

        for (label i = 0; i < nNodes_[dimj]; i++)
        {
            pos[dimj] = i;

            if (dimj == 0)
            {
                labelList posM(nDims_,0);
                for (label mi = 0; mi < pos.size(); mi++)
                {
                    posM[mi] = pos[mi];
                }
                Yold(i,0) = moments_(posM);
            }
            else
            {
                Yold(i,0) =
                    conditionalMoments_[dimi][dimj - 1](pos);
            }
        }

        labelList posVR(max(1, dimj), 0);
        if (dimj != 0)
        {
            for (label i = 0; i < posVR.size(); i++)
            {
                posVR[i] = pos[i];
            }
        }

        scalarRectangularMatrix Ynew = invVR_[dimj](posVR)*Yold;

        for (label i = 0; i < nNodes_[dimj]; i++)
        {
            pos[dimj] = i;
            conditionalMoments_[dimi][dimj](pos) = Ynew(i,0);
        }
    }
}

void Foam::multivariateMomentInversion::setVR
(
    const label dimi,
    labelList& pos,
    label ai
)
{
    if (ai < dimi)
    {
        for (label i = 0; i < nNodes_[ai]; i++)
        {
            pos[ai] = i;

            setVR(dimi,pos,ai + 1);
        }
    }
    else
    {
        scalarDiagonalMatrix x(nNodes_[dimi], 0.0);
        scalarSquareMatrix invR(nNodes_[dimi], 0.0);

        for (label i = 0; i < nNodes_[dimi]; i++)
        {
            pos[dimi] = i;
            x[i] = abscissae_[dimi](pos);
            invR[i][i] = 1.0/weights_[dimi](pos);
        }

        Vandermonde V(x);

        scalarSquareMatrix invV(V.inv());

        labelList posVR(max(1, dimi), 0);

        if (dimi > 0)
        {
            for (label ai = 0; ai < pos.size(); ai++)
            {
                posVR[ai] = pos[ai];
            }
        }

        invVR_[dimi](posVR) = invR*invV;
    }
}

void Foam::multivariateMomentInversion::cycleAlphaWheeler
(
    const label dimi,
    label ai,
    labelList& pos
)
{
    if (ai < dimi)
    {
        for (label i = 0; i < nNodes_[ai]; i++)
        {
            pos[ai] = i;

            cycleAlphaWheeler(dimi,ai+1,pos);
        }

        return;
    }
    else
    {
        univariateMomentSet momentsToInvert
        (
            2*nNodes_[dimi],
            0.0,
            "Gauss",
            support_[dimi]
        );

        forAll(momentsToInvert, mi)
        {
            pos[dimi] = mi;
            momentsToInvert[mi] =
                conditionalMoments_[dimi][dimi - 1](pos);
        }

        momentsToInvert.invert();

        for (label nodei = 0; nodei < nNodes_[ai]; nodei++)
        {
            pos[dimi] = nodei;

            weights_[dimi](pos) =
                momentsToInvert.weights()[nodei];

            abscissae_[dimi](pos) =
                momentsToInvert.abscissae()[nodei];
        }
        return;
    }
}

void Foam::multivariateMomentInversion::reset()
{
    forAll(moments_, mi)
    {
        moments_[mi] = 0.0;
    }

    forAll(invVR_, dimi)
    {
        forAll(invVR_[dimi], ai)
        {
            for (label i = 0; i < nNodes_[dimi]; i++)
            {
                for (label j = 0; j < nNodes_[dimi]; j++)
                {
                    invVR_[dimi][ai](i,j) = 0.0;
                }
            }
        }

        forAll(conditionalMoments_[dimi], dimj)
        {
            forAll(conditionalMoments_[dimi][dimj], ai)
            {
                conditionalMoments_[dimi][dimj][ai] = 0.0;
            }
        }
    }

    forAll(abscissae_, dimi)
    {
        forAll(abscissae_[dimi], ai)
        {
            weights_[dimi][ai] = 0.0;
            abscissae_[dimi][ai] = 0.0;
        }
    }
}


// ************************************************************************* //
