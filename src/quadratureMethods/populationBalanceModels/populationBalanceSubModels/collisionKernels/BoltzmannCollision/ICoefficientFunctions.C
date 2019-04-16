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

// Zero order
IFuncHeader(0,0,0)
{
    Is(0)[celli] = 0.0;
}

// First order
IFuncHeader(0,0,1)
{
    Is(0,0,1)[celli] = -(omega/2.0)*g.z();
}

IFuncHeader(0,1,0)
{
    Is(0,1)[celli] = -(omega/2.0)*g.y();
}

IFuncHeader(1,0,0)
{
    Is(1)[celli] = -(omega/2.0)*g.x();
}

// Second order
IFuncHeader(0,0,2)
{
    Is(0,0,2)[celli] =
        (sqr(omega)/12.0)*sqr(magg)
      + (sqr(omega)/4.0)*sqr(g.z())
      - omega*g.z()*v.z();
}

IFuncHeader(0,1,1)
{
    Is(0,1,1)[celli] =
        (sqr(omega)/4.0)*g.y()*g.z()
      - (omega/2.0)*(v.y()*g.z() + g.y()*v.z());
}

IFuncHeader(1,0,1)
{
    Is(1,0,1)[celli] =
        (sqr(omega)/4.0)*g.x()*g.z()
      - (omega/2.0)*(v1*g.z() + g.x()*v.z());
}

IFuncHeader(1,1,0)
{
    Is(1,1)[celli] =
        (sqr(omega)/4.0)*g.x()*g.y()
      - (omega/2.0)*(v.x()*g.y() + g.x()*v.y());
}

IFuncHeader(0,2,0)
{
    Is(0,2)[celli] =
        (sqr(omega)/12.0)*sqr(magg)
      + (sqr(omega)/4.0)*sqr(g.y())
      - omega*g.y()*v.y();
}

IFuncHeader(2,0,0)
{
    Is(2)[celli] =
        (sqr(omega)/12.0)*sqr(magg)
      + (sqr(omega)/4.0)*sqr(g.x())
      - omega*g.x()*v1;
}

// Third order
IFuncHeader(0,0,3)
{
    Is(0,0,3)[celli] =
      - (pow3(omega)/8.0)*(sqr(magg) + sqr(g.z()))*g.z()
      + (sqr(omega)/4.0)*(sqr(magg) + 3.0*sqr(g.z()))*v.z()
      - (1.5*omega)*g.z()*sqr(v.z());
}

IFuncHeader(0,1,2)
{
    Is(0,1,2)[celli] =
      - pow3(omega)/8.0*(sqr(magg) + 3.0*sqr(g.z()))*g.y()
      + sqr(omega)/2.0*g.z()*g.y()*v.z()
      + sqr(omega)/12.0*(sqr(magg) + 3.0*sqr(g.z()))*v.y()
      - omega/2.0*g.y()*sqr(v.z())
      - omega*g.z()*v.z()*v.y();
}

IFuncHeader(0,2,1)
{
    Is(0,2,1)[celli] =
      - pow3(omega)/8.0*(sqr(magg) + 3.0*sqr(g.y()))*g.z()
      + sqr(omega)/2.0*g.y()*g.z()*v.y()
      + sqr(omega)/12.0*(sqr(magg) + 3.0*sqr(g.y()))*v.z()
      - omega/2.0*g.z()*sqr(v.y())
      - omega*g.y()*v.y()*v.z();
}

IFuncHeader(0,3,0)
{
    Is(0,3)[celli] =
      - (pow3(omega)/8.0)*(sqr(magg) + sqr(g.y()))*g.y()
      + (sqr(omega)/4.0)*(sqr(magg) + 3.0*sqr(g.y()))*v.y()
      - (1.5*omega)*g.y()*sqr(v.y());
}

IFuncHeader(1,0,2)
{
    Is(1,0,2)[celli] =
      - pow3(omega)/8.0*(sqr(magg) + 3.0*sqr(g.z()))*g.x()
      + sqr(omega)/2.0*g.z()*g.x()*v.z()
      + sqr(omega)/12.0*(sqr(magg) + 3.0*sqr(g.z()))*v.x()
      - omega/2.0*g.x()*sqr(v.z())
      - omega*g.z()*v.z()*v.x();
}

IFuncHeader(1,1,1)
{
    Is(1,1,1)[celli] =
      - pow3(omega)/8.0*g.x()*g.y()*g.z()
      + sqr(omega)/4.0
       *(
            g.x()*g.y()*v.z()
          + g.y()*g.z()*v.x()
          + g.z()*g.x()*v.y()
        )
      - omega/2.0
       *(
            g.x()*v.y()*v.z()
          + g.y()*v.z()*v.x()
          + g.z()*v.x()*v.y()
        );
}

IFuncHeader(1,2,0)
{
    Is(1,2)[celli] =
      - pow3(omega)/8.0*(sqr(magg) + 3.0*sqr(g.y()))*g.x()
      + sqr(omega)/2.0*g.y()*g.x()*v.y()
      + sqr(omega)/12.0*(sqr(magg) + 3.0*sqr(g.y()))*v.x()
      - omega/2.0*g.x()*sqr(v.y())
      - omega*g.y()*v.y()*v.x();
}

IFuncHeader(2,0,1)
{
    Is(2,0,1)[celli] =
      - pow3(omega)/8.0*(sqr(magg) + 3.0*sqr(g.x()))*g.z()
      + sqr(omega)/2.0*g.x()*g.z()*v.x()
      + sqr(omega)/12.0*(sqr(magg) + 3.0*sqr(g.x()))*v.z()
      - omega/2.0*g.z()*sqr(v.x())
      - omega*g.x()*v.x()*v.z();
}

IFuncHeader(2,1,0)
{
    Is(2,1)[celli] =
      - pow3(omega)/8.0*(sqr(magg) + 3.0*sqr(g.x()))*g.y()
      + sqr(omega)/2.0*g.x()*g.y()*v.x()
      + sqr(omega)/12.0*(sqr(magg) + 3.0*sqr(g.x()))*v.y()
      - omega/2.0*g.y()*sqr(v.x())
      - omega*g.x()*v.x()*v.y();
}

IFuncHeader(3,0,0)
{
    Is(3)[celli] =
      - (pow3(omega)/8.0)*(sqr(magg) + sqr(g.x()))*g.x()
      + (sqr(omega)/4.0)*(sqr(magg) + 3.0*sqr(g.x()))*v1
      - (1.5*omega)*g.x()*sqr(v1);
}

// Fourth order
IFuncHeader(0,0,4)
{
    Is(0,0,4)[celli] =
        (pow4(omega)/80.0)*(pow4(magg) + 10.0*sqr(magg)*sqr(g.z()) + 5.0*pow4(g.z()))
      - (pow3(omega)/2.0)*(sqr(magg) + sqr(g.z()))*g.z()*v.z()
      + (sqr(omega)/2.0)*(sqr(magg) + 3.0*sqr(g.z()))*sqr(v.z())
      - 2.0*omega*g.z()*pow3(v.z());
}

// IFuncHeader(0,1,3)

// IFuncHeader(0,2,2)

// IFuncHeader(0,3,1)

IFuncHeader(0,4,0)
{
    Is(0,4)[celli] =
        (pow4(omega)/80.0)*(pow4(magg) + 10.0*sqr(magg)*sqr(g.y()) + 5.0*pow4(g.y()))
      - (pow3(omega)/2.0)*(sqr(magg) + sqr(g.y()))*g.y()*v.y()
      + (sqr(omega)/2.0)*(sqr(magg) + 3.0*sqr(g.y()))*sqr(v.y())
      - 2.0*omega*g.y()*pow3(v.y());
}

// IFuncHeader(1,0,3)

// IFuncHeader(2,0,2)

// IFuncHeader(2,2,0)

// IFuncHeader(1,3,0)

// IFuncHeader(3,0,1)

// IFuncHeader(3,1,0)


IFuncHeader(4,0,0)
{
    Is(4)[celli] =
        (pow4(omega)/80.0)*(pow4(magg) + 10.0*sqr(magg)*sqr(g.x()) + 5.0*pow4(g.x()))
      - (pow3(omega)/2.0)*(sqr(magg) + sqr(g.x()))*g.x()*v1
      + (sqr(omega)/2.0)*(sqr(magg) + 3.0*sqr(g.x()))*sqr(v1)
      - 2.0*omega*g.x()*pow3(v1);
}


// // Fith order
// IFuncHeader(0,0,5)

// IFuncHeader(0,1,4)

// IFuncHeader(0,2,3)

// IFuncHeader(0,3,2)

// IFuncHeader(0,4,1)

// IFuncHeader(0,5,0)

// IFuncHeader(1,0,4)

// IFuncHeader(1,4,0)

// IFuncHeader(2,0,3)

// IFuncHeader(2,3,0)

// IFuncHeader(3,0,2)

// IFuncHeader(3,2,0)

// IFuncHeader(4,0,1)

// IFuncHeader(4,1,0)

// IFuncHeader(5,0,0)


// // Sixth order

// IFuncHeader(0,1,5)

// IFuncHeader(0,2,4)

// IFuncHeader(0,4,2)

// IFuncHeader(0,5,1)

// IFuncHeader(1,0,5)

// IFuncHeader(1,5,0)

// IFuncHeader(2,0,4)

// IFuncHeader(2,4,0)

// IFuncHeader(4,0,2)

// IFuncHeader(4,2,0)

// IFuncHeader(5,0,1)

// IFuncHeader(5,1,0)

// IFuncHeader(0,2,5)

// IFuncHeader(0,5,2)

// IFuncHeader(2,0,5)

// IFuncHeader(2,5,0)

// IFuncHeader(5,0,2)

// IFuncHeader(5,2,0)

// ************************************************************************* //
