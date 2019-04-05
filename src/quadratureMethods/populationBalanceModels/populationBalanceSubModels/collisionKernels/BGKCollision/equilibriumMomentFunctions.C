/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Alberto Passalacqua
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
#include "BGKCollision.H"

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment000
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0)[celli] = m0;
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment100
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(1)[celli] = m0*u;
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment010
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,1)[celli] = m0*v;
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment001
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,0,1)[celli] = m0*w;
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment200
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(2)[celli] = m0*(sigma.xx() + sqr(u));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment020
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,2)[celli] = m0*(sigma.yy() + sqr(v));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment002
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,0,2)[celli] = m0*(sigma.zz() + sqr(w));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment300
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(3)[celli] = m0*(3.0*sigma.xx()*u + pow3(u));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment030
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,3)[celli] = m0*(3.0*sigma.yy()*v + pow3(v));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment003
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,0,3)[celli] = m0*(3.0*sigma.zz()*w + pow3(w));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment400
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(4)[celli] =
        m0*(6.0*sqr(u)*sigma.xx() + 3.0*sqr(sigma.xx()) + pow4(u));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment040
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,4)[celli] =
        m0*(6.0*sqr(v)*sigma.yy() + 3.0*sqr(sigma.yy()) + pow4(v));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment004
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,0,4)[celli] =
        m0*(6.0*sqr(w)*sigma.zz() + 3.0*sqr(sigma.zz()) + pow4(w));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment500
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(5)[celli] =
        m0*(15.0*u*sqr(sigma.xx()) + 10.0*sigma.xx()*pow3(u) + pow5(u));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment050
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,5)[celli] =
        m0*(15.0*v*sqr(sigma.yy()) + 10.0*sigma.yy()*pow3(v) + pow5(v));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment005
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,0,5)[celli] =
        m0*(15.0*w*sqr(sigma.zz()) + 10.0*sigma.zz()*pow3(w) + pow5(w));
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment110
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(1,1,0)[celli] = m0*(sigma.xy() + u*v);
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment101
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(1,0,1)[celli] = m0*(sigma.xz() + u*w);
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment011
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,1,1)[celli] = m0*(sigma.yz() + v*w);
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment210
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(2,1,0)[celli] = m0*(v*sqr(u) + 2.0*sigma.xy()*u + sigma.xx()*v);
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment201
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(2,0,1)[celli] = m0*(w*sqr(u) + 2.0*sigma.xz()*u + sigma.xx()*w);
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment021
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,2,1)[celli] = m0*(w*sqr(v) + 2.0*sigma.yz()*v + sigma.yy()*w);
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment120
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(1,2,0)[celli] = m0*(u*sqr(v) + 2.0*sigma.xy()*v + sigma.yy()*u);
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment102
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(1,0,2)[celli] = m0*(u*sqr(w) + 2.0*sigma.xz()*w + sigma.zz()*u);
}

void
Foam::populationBalanceSubModels::collisionKernels::BGKCollision::moment012
(
    mappedPtrList<volScalarField>& moments,
    const label celli,
    const scalar& m0,
    const scalar& u,
    const scalar& v,
    const scalar& w,
    const symmTensor& sigma
)
{
    moments(0,1,2)[celli] = m0*(v*sqr(w) + 2.0*sigma.yz()*w + sigma.zz()*v);
}
// ************************************************************************* //
