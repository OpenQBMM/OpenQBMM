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
IzFuncHeader(0,0,0)
{
    Is(0)[0][celli] = 0.0;
}

// First order
IzFuncHeader(0,0,1)
{
    Is(0,0,1)[0][celli] = (4.0*omega/15.0)*g.x()*g.z()
}

IzFuncHeader(0,1,0)
{
    Is(0,1)[0][celli] = (4.0*omega/15.0)*g.x()*g.y();
}

IzFuncHeader(1,0,0)
{
    Is(1)[0][celli] = (2.0*omega/15.0)*(sqr(magg) + 2.0*sqr(g.x()));
}

// Second order
IzFuncHeader(0,0,2)
{
    Is(0,0,2)[0][celli] =
      - (2.0*sqr(omega)/35.0)*(sqr(magg) + 2.0*sqr(g.z()))*g.x()
      + (8.0*omega/15.0)*g.x()*g.z()*v.z();
}

IzFuncHeader(0,1,1)
{
    Is(0,1,1)[0][celli] =
      - (4.0*sqr(omega)/35.0)*g.x()*g.y()*g.z()
      + (4.0*omega/15.0)*(g.x()*g.z()*v.y() + g.x()*g.y()*v.z());
}

IzFuncHeader(1,0,1)
{
    Is(1,0,1)[0][celli] =
      - (2.0*sqr(omega)/35.0)*(sqr(magg) + 2.0*sqr(g.z()))*(g.x() + v.z())
      + (4.0*omega/15.0)*g.x()*g.z()*v.x();
}

IzFuncHeader(1,1,0)
{
    Is(1,1)[0][celli] =
      - (2.0*sqr(omega)/35.0)*(sqr(magg) + 2.0*sqr(g.y()))*(g.x() + v.y())
      + (4.0*omega/15.0)*g.x()*g.y()*v.x();
}

IzFuncHeader(0,2,0)
{
    Is(0,2)[0][celli] =
      - (2.0*sqr(omega)/35.0)*(sqr(magg) + 2.0*sqr(g.y()))*g.x()
      + (8.0*omega/15.0)*g.x()*g.y()*v.y();
}

IzFuncHeader(2,0,0)
{
    Is(2)[0][celli] =
      - (2.0*sqr(omega)/35.0)*(3.0*sqr(magg) + 2.0*sqr(g.x()))*g.x()
      + (4.0*omega/15.0)*(sqr(magg) + 2.0*sqr(g.x()))*v.x();
}

// Third order
IzFuncHeader(0,0,3)
{
    Is(0,0,3)[0][celli] =
}

IzFuncHeader(0,1,2)
{
    Is(0,1,2)[0][celli] =
}

IzFuncHeader(0,2,1)
{
    Is(0,2,1)[0][celli] =
}

IzFuncHeader(0,3,0)
{
    Is(0,3)[0][celli] =
}

IzFuncHeader(1,0,2)
{
    Is(1,0,2)[0][celli] =
}

IzFuncHeader(1,1,1)
{
    Is(1,1,1)[0][celli] =
}

IzFuncHeader(1,2,0)
{
    Is(1,2)[0][celli] =
}

IzFuncHeader(2,0,1)
{
    Is(2,0,1)[0][celli] =
}

IzFuncHeader(2,1,0)
{
    Is(2,1)[0][celli] =
}

IzFuncHeader(3,0,0)
{
    Is(3)[0][celli] =
}

// Fourth order
IzFuncHeader(0,0,4)
{
    Is(0,0,4)[celli] =
}

// IzFuncHeader(0,1,3)

// IzFuncHeader(0,2,2)

// IzFuncHeader(0,3,1)

IzFuncHeader(0,4,0)
{
    Is(0,4)[0][celli] =
}

// IzFuncHeader(1,0,3)

// IzFuncHeader(2,0,2)

// IzFuncHeader(2,2,0)

// IzFuncHeader(1,3,0)

// IzFuncHeader(3,0,1)

// IzFuncHeader(3,1,0)


IzFuncHeader(4,0,0)
{
    Is(4)[0][celli] =
        (pow4(omega)/80.0)*(pow4(magg) + 10.0*sqr(magg)*sqr(g.x()) + 5.0*pow4(g.x()))
      - (pow3(omega)/2.0)*(sqr(magg) + sqr(g.x()))*g.x()*v1
      + (sqr(omega)/2.0)*(sqr(magg) + 3.0*sqr(g.x()))*sqr(v1)
      - 2.0*omega*g.x()*pow3(v1);
}


// // Fith order
// IzFuncHeader(0,0,5)

// IzFuncHeader(0,1,4)

// IzFuncHeader(0,2,3)

// IzFuncHeader(0,3,2)

// IzFuncHeader(0,4,1)

// IzFuncHeader(0,5,0)

// IzFuncHeader(1,0,4)

// IzFuncHeader(1,4,0)

// IzFuncHeader(2,0,3)

// IzFuncHeader(2,3,0)

// IzFuncHeader(3,0,2)

// IzFuncHeader(3,2,0)

// IzFuncHeader(4,0,1)

// IzFuncHeader(4,1,0)

// IzFuncHeader(5,0,0)


// // Sixth order

// IzFuncHeader(0,1,5)

// IzFuncHeader(0,2,4)

// IzFuncHeader(0,4,2)

// IzFuncHeader(0,5,1)

// IzFuncHeader(1,0,5)

// IzFuncHeader(1,5,0)

// IzFuncHeader(2,0,4)

// IzFuncHeader(2,4,0)

// IzFuncHeader(4,0,2)

// IzFuncHeader(4,2,0)

// IzFuncHeader(5,0,1)

// IzFuncHeader(5,1,0)

// IzFuncHeader(0,2,5)

// IzFuncHeader(0,5,2)

// IzFuncHeader(2,0,5)

// IzFuncHeader(2,5,0)

// IzFuncHeader(5,0,2)

// IzFuncHeader(5,2,0)

// ************************************************************************* //
