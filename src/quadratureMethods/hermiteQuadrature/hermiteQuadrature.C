/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2010-2014 by Bo Kong and 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2019-2020 Alberto Passalacqua
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

\*----------------------------------------------------------------------------*/

#include "hermiteQuadrature.H"
#include "scalar.H"
#include "scalarMatrices.H"
#include "EigenMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

const Foam::scalar Foam::hermiteQuadrature::thetaLimit_ = 1.0E-7;

Foam::hermiteQuadrature::hermiteQuadrature
(
    const label& nDim,
    const label& nOrder
)
:
    nDim_(nDim),
    nOrder_(nOrder),
    nTolNodes_(pow(nOrder_, nDim)),
    origWei_(nTolNodes_, Zero),
    origAbs_(nTolNodes_, Zero),
    resAbs_(nTolNodes_, Zero)
{

    if((nOrder_ <= 0))
    {
        FatalErrorIn("Foam::hermiteQuadrature\n" )
            << "parameter(s) out of range ! "
            << abort(FatalError);
    }

    scalarRectangularMatrix ab(nOrder_, 2, scalar(0));

    for(label i = 0; i < nOrder_; i++)
    {
        ab[i][1]= scalar(i);
    }

    scalarSquareMatrix z(nOrder_, Zero);

    for (label i = 0; i < nOrder_ - 1; i++)
    {
        z[i][i] = ab[i][0];
        z[i][i+1] = Foam::sqrt(ab[i+1][1]);
        z[i+1][i] = z[i][i+1];
    }

    z[nOrder_-1][nOrder_-1] = ab[nOrder_-1][0];

    EigenMatrix<scalar> zEig(z);

    scalarList herWei_(nOrder_, Zero);
    scalarList herAbs_(nOrder_, Zero);

    forAll(herWei_,i)
    {
        herWei_[i] = sqr(zEig.EVecs()[0][i]);
        herAbs_[i] = zEig.EValsRe()[i];
    }

    scalar wtotal = sum(herWei_) ;

    forAll(herWei_,i)
    {
        herWei_[i] = herWei_[i]/wtotal;
    }

    if (nDim_ == 1)
    {
        forAll(origWei_,i)
        {
            origWei_[i] = herWei_[i];
            origAbs_[i] = vector(herAbs_[i], 0, 0);
        }
    }

    maxAbs_ = max(herAbs_);

    if (nDim_ == 2)
    {
        for(label i = 0; i < nOrder_; i++)
        {
            for(label j = 0; j < nOrder_; j++)
            {
                label p = i*nOrder_ + j;
                origWei_[p] = herWei_[i]*herWei_[j];
                origAbs_[p] = vector(herAbs_[i], herAbs_[j], 0);

            }
        }
    }

    if (nDim_ == 3)
    {
        for(label i = 0; i < nOrder_; i++)
        {
            for(label j = 0; j < nOrder_; j++)
            {
                for(label k = 0; k < nOrder_; k++)
                {
                    label p = i*nOrder_*nOrder_ + j*nOrder_ + k;
                    origWei_[p] = herWei_[i]*herWei_[j]*herWei_[k];
                    origAbs_[p] = vector(herAbs_[i], herAbs_[j], herAbs_[k]);
                }
            }
        }
    }

    return ;
}

void Foam::hermiteQuadrature::calcHermiteQuadrature
(
    const vector& mu,
    const symmTensor& Pp
)
{
    if (tr(Pp) > thetaLimit_)
    {
        tensor spM(tensor::zero);

        if( nDim_ == 3)
        {
            scalarSquareMatrix z(3, Zero);

            z[0][0] = Pp.xx();
            z[0][1] = Pp.xy();
            z[0][2] = Pp.xz();
            z[1][0] = Pp.xy();
            z[1][1] = Pp.yy();
            z[1][2] = Pp.yz();
            z[2][0] = Pp.xz();
            z[2][1] = Pp.yz();
            z[2][2] = Pp.zz();

            EigenMatrix<scalar> zEig(z);
            const scalarDiagonalMatrix& e(zEig.EValsRe());
            const scalarSquareMatrix& ev(zEig.EVecs());

            scalarSquareMatrix E(3, Zero);

            forAll(e,i)
            {
                if(e[i] >= 0)
                {
                    E[i][i] = sqrt(e[i]);
                }
                else
                {
                    E[i][i] = 0.0;
                }
            }

            z = Zero;
            multiply(z, ev, E);
            forAll(spM,i) spM[i] = z[label(i/3)][i % 3];

        }
        else
        {
            if(nDim_==2)
            {
                scalarSquareMatrix z(2, Zero);
                z[0][0] = Pp.xx();
                z[0][1] = Pp.xy();
                z[1][0] = Pp.xy();
                z[1][1] = Pp.yy();

                EigenMatrix<scalar> zEig(z);
                const scalarDiagonalMatrix& e(zEig.EValsRe());
                const scalarSquareMatrix& ev(zEig.EVecs());

                scalarSquareMatrix E(2,Zero);

                forAll(e, i)
                {
                    if(e[i] >=0)
                    {
                        E[i][i] = sqrt(e[i]);
                    }
                    else
                    {
                        E[i][i] = 0.0;
                    }
                }

                z = Zero;
                multiply(z, ev, E);

                spM.xx() = z[0][0];
                spM.xy() = z[0][1];
                spM.yx() = z[1][0];
                spM.yy() = z[1][1];

            }
            else
            {
                spM.xx() = sqrt(Pp.xx());
            }
        }

        resAbs_ = (spM & origAbs_) + mu;

    }
    else
    {
        resAbs_ = mu ;
    }

    return ;
}

// ************************************************************************* //
