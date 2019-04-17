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
#include "BoltzmannCollision.H"

// Zero order
IzFuncHeader(0,0,0)
{
    Iz(0) = 0.0;
}

// First order
IzFuncHeader(0,0,1)
{
    Iz(0,0,1) = -(2.0*omegaPow[1]/15.0)*(gMagPow[2] + 2.0*gPow[2].z());
}

IzFuncHeader(0,1,0)
{
    Iz(0,1) = (4.0*omegaPow[1]/15.0)*gPow[1].y()*gPow[1].z();
}

IzFuncHeader(1,0,0)
{
    Iz(1) = (4.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].z();
}

// Second order
IzFuncHeader(0,0,2)
{
    Iz(0,0,2) =
      - (2.0*omegaPow[2]/35.0)*(3.0*gMagPow[2] + 2.0*gPow[2].z())*gPow[1].z()
      + (4.0*omegaPow[1]/15.0)*(gMagPow[2] + 2.0*gPow[2].z())*vPow[1].z();
}

IzFuncHeader(0,1,1)
{
    Iz(0,1,1) =
      - (2.0*omegaPow[2]/35.0)*(gMagPow[2] + 2.0*gPow[2].z())*gPow[1].y()
      + (4.0*omegaPow[1]/15.0)*gPow[1].y()*gPow[1].z()*vPow[1].z()
      + (2.0*omegaPow[1]/15.0)*(gMagPow[2] + 2.0*gPow[2].z())*vPow[1].y();
}

IzFuncHeader(1,0,1)
{
    Iz(1,0,1) =
      - (2.0*omegaPow[2]/35.0)*(gMagPow[2] + 2.0*gPow[2].z())*gPow[1].x()
      + (4.0*omegaPow[1]/15.0)*gPow[1].z()*gPow[1].x()*vPow[1].z()
      + (2.0*omegaPow[1]/15.0)*(gMagPow[2] + 2.0*gPow[2].z())*vPow[1].x();
}

IzFuncHeader(1,1,0)
{
    Iz(1,1) =
      - (4.0*omegaPow[2]/35.0)*gPow[1].x()*gPow[1].y()*gPow[1].z()
      + (4.0*omegaPow[1]/15.0)
       *gPow[1].z()*(gPow[1].x()*vPow[1].y() + gPow[1].y()*vPow[1].x());
}

IzFuncHeader(0,2,0)
{
    Iz(0,2) =
      - (2.0*omegaPow[2]/35.0)*(gMagPow[2] + 2.0*gPow[2].y())*gPow[1].z()
      + (8.0*omegaPow[1]/15.0)*gPow[1].y()*gPow[1].z()*vPow[1].y();
}

IzFuncHeader(2,0,0)
{
    Iz(2) =
      - (2.0*omegaPow[2]/35.0)*(gMagPow[2] + 2.0*gPow[2].x())*gPow[1].z()
      + (8.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].z()*vPow[1].x();
}

// Third order
IzFuncHeader(0,0,3)
{
    Iz(0,0,3) =
        (2.0*omegaPow[2]/315.0)
       *(3.0*gMagPow[4] + 24.0*gMagPow[2]*gPow[2].z() + 8.0*gPow[4].z())
      - (6.0*omegaPow[2]/35.0)
       *(3.0*gMagPow[2] + 2.0*gPow[2].z())*gPow[1].z()*vPow[1].z()
      + (2.0*omegaPow[1]/5.0)*(gMagPow[2] + 2.0*gPow[2].z())*vPow[2].z();
}

// IzFuncHeader(0,1,2)

// IzFuncHeader(0,2,1)


IzFuncHeader(0,3,0)
{
    Iz(0,3) =
        (8.0*omegaPow[3]/315.0)
       *(3.0*gMagPow[2] + 2.0*gPow[2].y())*gPow[1].z()*gPow[1].y()
      - (6.0*omegaPow[2]/35.0)
       *(gMagPow[2] + 2.0*gPow[2].y())*gPow[1].z()*vPow[1].y()
      + (4.0*omegaPow[1]/5.0)*gPow[1].z()*gPow[1].y()*vPow[2].y();
}

// IzFuncHeader(1,0,2)

// IzFuncHeader(1,1,1)

// IzFuncHeader(1,2,0)

// IzFuncHeader(2,0,1)

// IzFuncHeader(2,1,0)
// {
//     Iz(2,1) =
//         (8.0*omegaPow[3]/315.0)*(3.0*gMagPow[2] + 2.0*gPow[2].x())*gPow[1].x()*gPow[1].y()
//       - (2.0*omegaPow[2]/35.0)
//        *(
//            2.0*(gMagPow[2] + 2.0*gPow[2].x())*gPow[1].y()*vPow[1].x()
//          + (3.0*gMagPow[2] + 2.0*gPow[2].x())*gPow[1].x()*vPow[1].y()
//         )
//       + (4.0*omegaPow[1]/15.0)
//        *(
//             (gMagPow[2] + 2.0*gPow[2].x())*vPow[1].x()*vPow[1].y()
//           + gPow[1].x()*gPow[1].y()*vPow[2].x()
//         );
// }


IzFuncHeader(3,0,0)
{
    Iz(3) =
        (8.0*omegaPow[3]/315.0)
       *(3.0*gMagPow[2] + 2.0*gPow[2].x())*gPow[1].z()*gPow[1].x()
      - (6.0*omegaPow[2]/35.0)
       *(gMagPow[2] + 2.0*gPow[2].x())*gPow[1].z()*vPow[1].x()
      + (4.0*omegaPow[1]/5.0)*gPow[1].x()*gPow[1].z()*vPow[2].x();
}

// Fourth order
IzFuncHeader(0,0,4)
{
    Iz(0,0,4) =
      - (2.0*omegaPow[4]/693.0)
       *(15.0*gMagPow[4] + 40.0*gMagPow[2]*gPow[2].z() + 8.0*gPow[4].z())
       *gPow[1].z()
      + (8.0*omegaPow[3]/315.0)
       *(3.0*gMagPow[4] + 24.0*gMagPow[2]*gPow[2].z() + 8.0*gPow[4].z())
       *vPow[1].z()
      - (12.0*omegaPow[2]/35.0)
       *(3.0*gMagPow[2] + 2.0*gPow[2].z())*gPow[1].z()*vPow[2].z()
      + (8.0*omegaPow[1]/15.0)*(gMagPow[2] + 2.0*gPow[2].z())*pow3(vPow[1].z());
}

IzFuncHeader(0,4,0)
{
    Iz(0,4) =
      - (2.0*omegaPow[4]/693.0)
       *(3.0*gMagPow[4] + 24.0*gMagPow[2]*gPow[2].y() + 8.0*gPow[4].y())*gPow[1].z()
      + (32.0*omegaPow[3]/315.0)
       *(3.0*gMagPow[2] + 2.0*gPow[2].y())*gPow[1].z()*gPow[1].y()*vPow[1].y()
      - (12.0*omegaPow[2]/35.0)
       *(gMagPow[2] + 2.0*gPow[2].y())*gPow[1].z()*vPow[2].y()
      + (16.0*omegaPow[1]/15.0)*gPow[1].z()*gPow[1].y()*vPow[3].y();
}

IzFuncHeader(4,0,0)
{
    Iz(4) =
      - (2.0*omegaPow[4]/693.0)
       *(3.0*gMagPow[4] + 24.0*gMagPow[2]*gPow[2].x() + 8.0*gPow[4].x())*gPow[1].z()
      + (32.0*omegaPow[3]/315.0)
       *(3.0*gMagPow[2] + 2.0*gPow[2].x())*gPow[1].z()*gPow[1].x()*vPow[1].x()
      - (12.0*omegaPow[2]/35.0)
       *(gMagPow[2] + 2.0*gPow[2].x())*gPow[1].z()*vPow[2].x()
      + (16.0*omegaPow[1]/15.0)*gPow[1].z()*gPow[1].x()*vPow[3].x();
}

// ************************************************************************* //
