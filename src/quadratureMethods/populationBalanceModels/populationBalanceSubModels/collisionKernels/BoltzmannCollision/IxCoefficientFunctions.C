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
IxFuncHeader(0,0,0)
{
    Ix(0) = 0.0;
}

// First order
IxFuncHeader(0,0,1)
{
    Ix(0,0,1) = (4.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].z();
}

IxFuncHeader(0,1,0)
{
    Ix(0,1) = (4.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].y();
}

IxFuncHeader(1,0,0)
{
    Ix(1) = (2.0*omegaPow[1]/15.0)*(gMagSqr + 2.0*gPow[2].x());
}

// Second order
IxFuncHeader(0,0,2)
{
    Ix(0,0,2) =
      - (2.0*omegaPow[2]/35.0)*(gMagSqr + 2.0*gPow[2].z())*gPow[1].x()
      + (8.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].z()*vPow[1].z();
}

IxFuncHeader(0,1,1)
{
    Ix(0,1,1) =
      - (4.0*omegaPow[2]/35.0)*gPow[1].x()*gPow[1].y()*gPow[1].z()
      + (4.0*omegaPow[1]/15.0)*(gPow[1].x()*gPow[1].z()*vPow[1].y() 
      + gPow[1].x()*gPow[1].y()*vPow[1].z());
}

IxFuncHeader(1,0,1)
{
    Ix(1,0,1) =
      - (2.0*omegaPow[2]/35.0)*(gMagSqr + 2.0*gPow[2].z())
        *(gPow[1].x() + vPow[1].z())
      + (4.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].z()*vPow[1].x();
}

IxFuncHeader(1,1,0)
{
    Ix(1,1) =
      - (2.0*omegaPow[2]/35.0)*(gMagSqr + 2.0*gPow[2].y())
        *(gPow[1].x() + vPow[1].y())
      + (4.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].y()*vPow[1].x();
}

IxFuncHeader(0,2,0)
{
    Ix(0,2) =
      - (2.0*omegaPow[2]/35.0)*(gMagSqr + 2.0*gPow[2].y())*gPow[1].x()
      + (8.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].y()*vPow[1].y();
}

IxFuncHeader(2,0,0)
{
    Ix(2) =
      - (2.0*omegaPow[2]/35.0)*(3.0*gMagSqr + 2.0*gPow[2].x())*gPow[1].x()
      + (4.0*omegaPow[1]/15.0)*(gMagSqr + 2.0*gPow[2].x())*vPow[1].x();
}

// Third order
IxFuncHeader(0,0,3)
{
    Ix(0,0,3) =
        (8.0*omegaPow[3]/315.0)
       *(3.0*gMagSqr + 2.0*gPow[2].z())*gPow[1].x()*gPow[1].z()
      - (6.0*omegaPow[2]/35.0)
       *(gMagSqr + 2.0*gPow[2].z())*gPow[1].x()*vPow[1].z()
      + (4.0*omegaPow[1]/5.0)*gPow[1].x()*gPow[1].z()*vPow[2].z();
}

IxFuncHeader(0,1,2)
{
    Ix(0,1,2) =
        (8.0*omegaPow[3]/315.0)
       *(3.0*gMagSqr + 2.0*gPow[2].z())*gPow[1].x()*gPow[1].y()
      - (2.0*omegaPow[2]/35.0)
       *(
           (gMagSqr + 2.0*gPow[2].z())*gPow[1].x()*vPow[1].y()
         + 4.0*gPow[1].x()*gPow[1].y()*gPow[1].z()*vPow[1].z()
        )
      + (4.0*omegaPow[1]/15.0)*gPow[1].x()*vPow[1].z()
       *(gPow[1].y()*vPow[1].z() + 2.0*gPow[1].z()*vPow[1].y());
}

IxFuncHeader(0,2,1)
{
    Ix(0,1,2) =
        (8.0*omegaPow[3]/315.0)
       *(3.0*gMagSqr + 2.0*gPow[2].y())*gPow[1].x()*gPow[1].z()
      - (2.0*omegaPow[2]/35.0)
       *(
           (gMagSqr + 2.0*gPow[2].y())*gPow[1].x()*vPow[1].z()
         + 4.0*gPow[1].x()*gPow[1].z()*gPow[1].y()*vPow[1].y()
        )
      + (4.0*omegaPow[1]/15.0)*gPow[1].x()*vPow[1].y()
       *(gPow[1].z()*vPow[1].y() + 2.0*gPow[1].y()*vPow[1].z());
}

IxFuncHeader(0,3,0)
{
    Ix(0,3) =
        (8.0*omegaPow[3]/315.0)
       *(3.0*gMagSqr + 2.0*gPow[2].y())*gPow[1].x()*gPow[1].y()
      - (6.0*omegaPow[2]/35.0)
       *(gMagSqr + 2.0*gPow[2].y())*gPow[1].x()*vPow[1].y()
      + (4.0*omegaPow[1]/5.0)*gPow[1].x()*gPow[1].y()*sqr(vPow[1].y());
}

IxFuncHeader(1,0,2)
{
    Ix(1,0,2) =
        (2.0*omegaPow[3]/315.0)
       *(
            sqr(gMagSqr)
          + 4.0*gMagSqr*(gPow[2].x() + gPow[2].z())
          + 8.0*gPow[2].x()*gPow[2].z()
        )
      - (2.0*omegaPow[2]/35.0)
       *(
           (gMagSqr + 2.0*gPow[2].z())*gPow[1].x()*vPow[1].x()
         + (gMagSqr + 2.0*gPow[2].x())*gPow[1].z()*vPow[1].z()
        )
      + (2.0*omegaPow[1]/15.0)
       *(
            (gMagSqr + 2.0*gPow[2].x())*vPow[2].z()
          + 4.0*gPow[1].x()*gPow[1].z()*vPow[1].x()*vPow[1].z()
        );
}

IxFuncHeader(1,1,1)
{
    Ix(1,1,1) =
        (8.0*omegaPow[3]/315.0)
       *(3.0*gMagSqr + 2.0*gPow[2].x())*gPow[1].y()*gPow[1].z()
      - (2.0*omegaPow[2]/35.0)
       *(
            2.0*gPow[1].x()*gPow[1].y()*gPow[1].z()*vPow[1].x()
          + (gMagSqr + 2.0*gPow[2].x())
           *(gPow[1].y()*vPow[1].z() + gPow[1].z()*vPow[1].y())
        )
      + (4.0*omegaPow[1]/15.0)
       *(
           2.0*gPow[1].x()*vPow[1].x()
          *(gPow[1].y()*vPow[1].z() + gPow[1].z()*vPow[1].y())
         + (gMagSqr + 2.0*gPow[2].x())*vPow[1].y()*vPow[1].z()
        );
}

IxFuncHeader(1,2,0)
{
    Ix(1,2,0) =
        (2.0*omegaPow[3]/315.0)
       *(
            sqr(gMagSqr)
          + 4.0*gMagSqr*(gPow[2].x() + gPow[2].y())
          + 8.0*gPow[2].x()*gPow[2].y()
        )
      - (2.0*omegaPow[2]/35.0)
       *(
           (gMagSqr + 2.0*gPow[2].y())*gPow[1].x()*vPow[1].x()
         + (gMagSqr + 2.0*gPow[2].x())*gPow[1].y()*vPow[1].y()
        )
      + (2.0*omegaPow[1]/15.0)
       *(
            (gMagSqr + 2.0*gPow[2].x())*vPow[2].y()
          + 4.0*gPow[1].x()*gPow[1].y()*vPow[1].x()*vPow[1].y()
        );
}

IxFuncHeader(2,0,1)
{
    Ix(2,0,1) =
        (8.0*omegaPow[3]/315.0)*(3.0*gMagSqr 
            + 2.0*gPow[2].x())*gPow[1].x()*gPow[1].z()
      - (2.0*omegaPow[2]/35.0)
       *(
           2.0*(gMagSqr + 2.0*gPow[2].x())*gPow[1].z()*vPow[1].x()
         + (3.0*gMagSqr + 2.0*gPow[2].x())*gPow[1].x()*vPow[1].z()
        )
      + (4.0*omegaPow[1]/15.0)
       *(
            (gMagSqr + 2.0*gPow[2].x())*vPow[1].x()*vPow[1].z()
          + gPow[1].x()*gPow[1].z()*vPow[2].x()
        );
}

IxFuncHeader(2,1,0)
{
    Ix(2,1) =
        (8.0*omegaPow[3]/315.0)*(3.0*gMagSqr + 2.0*gPow[2].x())*gPow[1].x()*gPow[1].y()
      - (2.0*omegaPow[2]/35.0)
       *(
           2.0*(gMagSqr + 2.0*gPow[2].x())*gPow[1].y()*vPow[1].x()
         + (3.0*gMagSqr + 2.0*gPow[2].x())*gPow[1].x()*vPow[1].y()
        )
      + (4.0*omegaPow[1]/15.0)
       *(
            (gMagSqr + 2.0*gPow[2].x())*vPow[1].x()*vPow[1].y()
          + gPow[1].x()*gPow[1].y()*vPow[2].x()
        );
}

IxFuncHeader(3,0,0)
{
    Ix(3) =
        (2.0*omegaPow[2]/315.0)
       *(3.0*sqr(gMagSqr) + 24.0*gMagSqr*gPow[2].x() + 8.0*gPow[4].x())
      - (6.0*omegaPow[2]/35.0)
       *(3.0*gMagSqr + 2.0*gPow[2].x())*gPow[1].x()*vPow[1].x()
      + (2.0*omegaPow[1]/5.0)*(gMagSqr + 2.0*gPow[2].x())*vPow[2].x();
}


// Fourth order
IxFuncHeader(0,0,4)
{
    Ix(0,0,4) =
      - (2.0*omegaPow[4]/693.0)
       *(3.0*sqr(gMagSqr) + 24.0*gMagSqr*gPow[2].z() + 8.0*gPow[4].z())*gPow[1].x()
      + (32.0*omegaPow[3]/315.0)
       *(3.0*gMagSqr + 2.0*gPow[2].z())*gPow[1].x()*gPow[1].z()*vPow[1].z()
      - (12.0*omegaPow[2]/35.0)
       *(gMagSqr + 2.0*gPow[2].z())*gPow[1].x()*vPow[2].z()
      + (16.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].z()*pow3(vPow[1].z());
}

IxFuncHeader(0,4,0)
{
    Ix(0,4) =
      - (2.0*omegaPow[4]/693.0)
       *(3.0*sqr(gMagSqr) + 24.0*gMagSqr*gPow[2].y() + 8.0*gPow[4].y())*gPow[1].x()
      + (32.0*omegaPow[3]/315.0)
       *(3.0*gMagSqr + 2.0*gPow[2].y())*gPow[1].x()*gPow[1].y()*vPow[1].y()
      - (12.0*omegaPow[2]/35.0)
       *(gMagSqr + 2.0*gPow[2].y())*gPow[1].x()*sqr(vPow[1].y())
      + (16.0*omegaPow[1]/15.0)*gPow[1].x()*gPow[1].y()*vPow[3].y();
}

IxFuncHeader(4,0,0)
{
    Ix(4) =
      - (2.0*omegaPow[4]/693.0)
       *(15.0*sqr(gMagSqr) + 40.0*gMagSqr*gPow[2].x() + 8.0*gPow[4].x())
       *gPow[1].x()
      + (8.0*omegaPow[3]/315.0)
       *(3.0*sqr(gMagSqr) + 24.0*gMagSqr*gPow[2].x() 
        + 8.0*gPow[4].x())*vPow[1].x()
      - (12.0*omegaPow[2]/35.0)
       *(3.0*gMagSqr + 2.0*gPow[2].x())*gPow[1].x()*vPow[2].x()
      + (8.0*omegaPow[1]/15.0)
       *(gMagSqr + 2.0*gPow[2].x())*vPow[3].x();
}

// ************************************************************************* //
