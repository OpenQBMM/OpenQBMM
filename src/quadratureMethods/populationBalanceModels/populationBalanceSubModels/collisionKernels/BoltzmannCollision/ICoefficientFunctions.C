/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Alberto Passalacqua
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
#include "BoltzmannCollision.H"

// Zero order
IFuncHeader(0,0,0)
{
    Is(0,0,0) = 0.0;
}

// First order
IFuncHeader(0,0,1)
{
    Is(0,0,1) = -(omegaPow[1]/2.0)*gPow[1].z();
}

IFuncHeader(0,1,0)
{
    Is(0,1,0) = -(omegaPow[1]/2.0)*gPow[1].y();
}

IFuncHeader(1,0,0)
{
    Is(1,0,0) = -(omegaPow[1]/2.0)*gPow[1].x();
}

// Second order
IFuncHeader(0,0,2)
{
    Is(0,0,2) =
        (omegaPow[2]/12.0)*(gMagSqr + 3.0*gPow[2].z())
      - omegaPow[1]*gPow[1].z()*vPow[1].z();
}

IFuncHeader(0,1,1)
{
    Is(0,1,1) =
        (omegaPow[2]/4.0)*gPow[1].y()*gPow[1].z()
      - (omegaPow[1]/2.0)*(vPow[1].y()*gPow[1].z() + gPow[1].y()*vPow[1].z());
}

IFuncHeader(1,0,1)
{
    Is(1,0,1) =
        (omegaPow[2]/4.0)*gPow[1].x()*gPow[1].z()
      - (omegaPow[1]/2.0)*(vPow[1].x()*gPow[1].z() + gPow[1].x()*vPow[1].z());
}

IFuncHeader(1,1,0)
{
    Is(1,1,0) =
        (omegaPow[2]/4.0)*gPow[1].x()*gPow[1].y()
      - (omegaPow[1]/2.0)*(vPow[1].x()*gPow[1].y() + gPow[1].x()*vPow[1].y());
}

IFuncHeader(0,2,0)
{
    Is(0,2,0) =
        (omegaPow[2]/12.0)*(gMagSqr + 3.0*gPow[2].y())
      - omegaPow[1]*gPow[1].y()*vPow[1].y();
}

IFuncHeader(2,0,0)
{
    Is(2,0,0) =
        (omegaPow[2]/12.0)*(gMagSqr + 3.0*gPow[2].x())
      - omegaPow[1]*gPow[1].x()*vPow[1].x();
}

// Third order
IFuncHeader(0,0,3)
{
    Is(0,0,3) =
      - (omegaPow[3]/8.0)*(gMagSqr + gPow[2].z())*gPow[1].z()
      + (omegaPow[2]/4.0)*(gMagSqr + 3.0*gPow[2].z())*vPow[1].z()
      - (1.5*omegaPow[1])*gPow[1].z()*vPow[2].z();
}

IFuncHeader(0,1,2)
{
    Is(0,1,2) =
      - omegaPow[3]/24.0*(gMagSqr + 3.0*gPow[2].z())*gPow[1].y()
      + omegaPow[2]/2.0*gPow[1].z()*gPow[1].y()*vPow[1].z()
      + omegaPow[2]/12.0*(gMagSqr + 3.0*gPow[2].z())*vPow[1].y()
      - omegaPow[1]/2.0*gPow[1].y()*vPow[2].z()
      - omegaPow[1]*gPow[1].z()*vPow[1].z()*vPow[1].y();
}

IFuncHeader(0,2,1)
{
    Is(0,2,1) =
      - omegaPow[3]/24.0*(gMagSqr + 3.0*gPow[2].y())*gPow[1].z()
      + omegaPow[2]/2.0*gPow[1].y()*gPow[1].z()*vPow[1].y()
      + omegaPow[2]/12.0*(gMagSqr + 3.0*gPow[2].y())*vPow[1].z()
      - omegaPow[1]/2.0*gPow[1].z()*vPow[2].y()
      - omegaPow[1]*gPow[1].y()*vPow[1].y()*vPow[1].z();
}

IFuncHeader(0,3,0)
{
    Is(0,3,0) =
      - (omegaPow[3]/8.0)*(gMagSqr + gPow[2].y())*gPow[1].y()
      + (omegaPow[2]/4.0)*(gMagSqr + 3.0*gPow[2].y())*vPow[1].y()
      - (1.5*omegaPow[1])*gPow[1].y()*vPow[2].y();
}

IFuncHeader(1,0,2)
{
    Is(1,0,2) =
      - omegaPow[3]/24.0*(gMagSqr + 3.0*gPow[2].z())*gPow[1].x()
      + omegaPow[2]/2.0*gPow[1].z()*gPow[1].x()*vPow[1].z()
      + omegaPow[2]/12.0*(gMagSqr + 3.0*gPow[2].z())*vPow[1].x()
      - omegaPow[1]/2.0*gPow[1].x()*vPow[2].z()
      - omegaPow[1]*gPow[1].z()*vPow[1].z()*vPow[1].x();
}

IFuncHeader(1,1,1)
{
    Is(1,1,1) =
      - omegaPow[3]/8.0*gPow[1].x()*gPow[1].y()*gPow[1].z()
      + omegaPow[2]/4.0
       *(
            gPow[1].x()*gPow[1].y()*vPow[1].z()
          + gPow[1].y()*gPow[1].z()*vPow[1].x()
          + gPow[1].z()*gPow[1].x()*vPow[1].y()
        )
      - omegaPow[1]/2.0
       *(
            gPow[1].x()*vPow[1].y()*vPow[1].z()
          + gPow[1].y()*vPow[1].z()*vPow[1].x()
          + gPow[1].z()*vPow[1].x()*vPow[1].y()
        );
}

IFuncHeader(1,2,0)
{
    Is(1,2,0) =
      - omegaPow[3]/24.0*(gMagSqr + 3.0*gPow[2].y())*gPow[1].x()
      + omegaPow[2]/2.0*gPow[1].y()*gPow[1].x()*vPow[1].y()
      + omegaPow[2]/12.0*(gMagSqr + 3.0*gPow[2].y())*vPow[1].x()
      - omegaPow[1]/2.0*gPow[1].x()*vPow[2].y()
      - omegaPow[1]*gPow[1].y()*vPow[1].y()*vPow[1].x();
}

IFuncHeader(2,0,1)
{
    Is(2,0,1) =
      - omegaPow[3]/24.0*(gMagSqr + 3.0*gPow[2].x())*gPow[1].z()
      + omegaPow[2]/2.0*gPow[1].x()*gPow[1].z()*vPow[1].x()
      + omegaPow[2]/12.0*(gMagSqr + 3.0*gPow[2].x())*vPow[1].z()
      - omegaPow[1]/2.0*gPow[1].z()*vPow[2].x()
      - omegaPow[1]*gPow[1].x()*vPow[1].x()*vPow[1].z();
}

IFuncHeader(2,1,0)
{
    Is(2,1,0) =
      - omegaPow[3]/24.0*(gMagSqr + 3.0*gPow[2].x())*gPow[1].y()
      + omegaPow[2]/2.0*gPow[1].x()*gPow[1].y()*vPow[1].x()
      + omegaPow[2]/12.0*(gMagSqr + 3.0*gPow[2].x())*vPow[1].y()
      - omegaPow[1]/2.0*gPow[1].y()*vPow[2].x()
      - omegaPow[1]*gPow[1].x()*vPow[1].x()*vPow[1].y();
}

IFuncHeader(3,0,0)
{
    Is(3,0,0) =
      - (omegaPow[3]/8.0)*(gMagSqr + gPow[2].x())*gPow[1].x()
      + (omegaPow[2]/4.0)*(gMagSqr + 3.0*gPow[2].x())*vPow[1].x()
      - (1.5*omegaPow[1])*gPow[1].x()*vPow[2].x();
}

// Fourth order
IFuncHeader(0,0,4)
{
    Is(0,0,4) =
        (omegaPow[4]/80.0)
       *((gPow[2] & gPow[2]) + 10.0*gMagSqr*gPow[2].z() + 5.0*gPow[4].z())
      - (omegaPow[3]/2.0)*(gMagSqr + gPow[2].z())*gPow[1].z()*vPow[1].z()
      + (omegaPow[2]/2.0)*(gMagSqr + 3.0*gPow[2].z())*vPow[2].z()
      - 2.0*omegaPow[1]*gPow[1].z()*vPow[3].z();
}

IFuncHeader(0,4,0)
{
    Is(0,4,0) =
        (omegaPow[4]/80.0)
       *((gPow[2] & gPow[2]) + 10.0*gMagSqr*gPow[2].y() + 5.0*gPow[4].y())
      - (omegaPow[3]/2.0)*(gMagSqr + gPow[2].y())*gPow[1].y()*vPow[1].y()
      + (omegaPow[2]/2.0)*(gMagSqr + 3.0*gPow[2].y())*vPow[2].y()
      - 2.0*omegaPow[1]*gPow[1].y()*vPow[3].y();
}

IFuncHeader(4,0,0)
{
    Is(4,0,0) =
        (omegaPow[4]/80.0)
       *((gPow[2] & gPow[2]) + 10.0*gMagSqr*gPow[2].x() + 5.0*gPow[4].x())
      - (omegaPow[3]/2.0)*(gMagSqr + gPow[2].x())*gPow[1].x()*vPow[1].x()
      + (omegaPow[2]/2.0)*(gMagSqr + 3.0*gPow[2].x())*vPow[2].x()
      - 2.0*omegaPow[1]*gPow[1].x()*vPow[3].x();
}

// ************************************************************************* //
