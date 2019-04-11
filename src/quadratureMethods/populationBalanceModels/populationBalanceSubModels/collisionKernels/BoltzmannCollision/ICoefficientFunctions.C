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
// {
//     Is(0,1,3)[celli] =
//         m0
//        *(
//             v*pow3(w)
//           + 3.0*sigma.yz()*sqr(w)
//           + 3.0*sigma.zz()*w*v
//           + 3.0*sigma.zz()*sigma.yz()
//         );
// }
//
IFuncHeader(0,2,2)
{
    Is(0,2,2)[celli] =
        m0
       *(
            2.0*sqr(sigma.yz())
          + 4.0*sigma.yz()*v*w
          + sqr(v*w)
          + sigma.zz()*sqr(v)
          + sigma.yy()*sqr(w)
          + sigma.yy()*sigma.zz()

        );
}
//
// IFuncHeader(0,3,1)
// {
//     Is(0,3,1)[celli] =
//         m0
//        *(
//             w*pow3(v)
//           + 3.0*sigma.yz()*sqr(v)
//           + 3.0*sigma.yy()*v*w
//           + 3.0*sigma.yy()*sigma.yz()
//         );
// }

IFuncHeader(0,4,0)
{
    Is(0,4)[celli] =
        (pow4(omega)/80.0)*(pow4(magg) + 10.0*sqr(magg)*sqr(g.y()) + 5.0*pow4(g.y()))
      - (pow3(omega)/2.0)*(sqr(magg) + sqr(g.y()))*g.y()*v.y()
      + (sqr(omega)/2.0)*(sqr(magg) + 3.0*sqr(g.y()))*sqr(v.y())
      - 2.0*omega*g.y()*pow3(v.y());
}

// IFuncHeader(1,0,3)
// {
//     Is(1,0,3)[celli] =
//         m0
//        *(
//             u*pow3(w)
//           + 3.0*sigma.xz()*sqr(w)
//           + 3.0*sigma.zz()*w*u
//           + 3.0*sigma.zz()*sigma.xz()
//         );
// }
//
// IFuncHeader(2,0,2)
// {
//     Is(2,0,2)[celli] =
//         m0
//        *(
//             2.0*sqr(sigma.xz())
//           + 4.0*sigma.xz()*u*w
//           + sqr(u*w)
//           + sigma.zz()*sqr(u) + sigma.xx()*sqr(w)
//           + sigma.xx()*sigma.zz()
//         );
// }
//
// IFuncHeader(2,2,0)
// {
//     Is(2,2)[celli] =
//         m0
//        *(
//             2.0*sqr(sigma.xy())
//           + 4.0*sigma.xy()*u*v
//           + sqr(u*v)
//           + sigma.yy()*sqr(u)
//           + sigma.xx()*sqr(v)
//           + sigma.xx()*sigma.yy()
//         );
// }
//
// IFuncHeader(1,3,0)
// {
//     Is(1,3)[celli] =
//         m0
//        *(
//             u*pow3(v)
//           + 3.0*sigma.xy()*sqr(v)
//           + 3.0*sigma.yy()*v*u
//           + 3.0*sigma.yy()*sigma.xy()
//         );
// }
//
// IFuncHeader(3,0,1)
// {
//     Is(3,0,1)[celli] =
//         m0
//        *(
//             w*pow3(u)
//           + 3.0*sigma.xz()*sqr(u)
//           + 3.0*sigma.xx()*u*w
//           + 3.0*sigma.xx()*sigma.xz()
//         );
// }
//
// IFuncHeader(3,1,0)
// {
//     Is(3,1)[celli] =
//         m0
//        *(
//             v*pow3(u)
//           + 3.0*sigma.xy()*sqr(u)
//           + 3.0*sigma.xx()*u*v
//           + 3.0*sigma.xx()*sigma.xy()
//         );
// }

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
// {
//     Is(0,0,5)[celli] =
//         m0*(15.0*w*sqr(sigma.zz()) + 10.0*sigma.zz()*pow3(w) + pow5(w));
// }
//
// IFuncHeader(0,1,4)
// {
//     Is(0,1,4)[celli] =
//         m0
//        *(
//             3.0*v*sqr(sigma.zz())
//           + 6.0*v*sigma.zz()*sqr(w)
//           + 12.0*sigma.yz()*sigma.zz()*w
//           + v*pow4(w)
//           + 4.0*sigma.yz()*pow3(w)
//         );
// }
//
// IFuncHeader(0,2,3)
// {
//     Is(0,2,3)[celli] =
//         m0
//        *(
//             sqr(v)*pow3(w)
//           + sigma.yy()*pow3(w)
//           + 6.0*sqr(sigma.yz())*w
//           + 6.0*sigma.yz()*v*sqr(w)
//           + 3.0*sigma.zz()*sqr(v)*w
//           + 6.0*sigma.yz()*sigma.zz()*v
//           + 3.0*sigma.yy()*sigma.zz()*w
//         );
// }
//
// IFuncHeader(0,3,2)
// {
//     Is(0,3,2)[celli] =
//         m0
//        *(
//             pow3(v)*sqr(w)
//           + 6.0*sqr(sigma.yz())*v
//           + sigma.zz()*pow3(v)
//           + 3.0*sigma.yy()*v*sqr(w)
//           + 6.0*sigma.yz()*sqr(v)*w
//           + 3.0*sigma.yy()*sigma.zz()*v
//           + 6.0*sigma.yy()*sigma.yz()*w
//         );
// }
//
// IFuncHeader(0,4,1)
// {
//     Is(0,4,1)[celli] =
//         m0
//        *(
//             3.0*w*sqr(sigma.yy())
//           + 6.0*w*sigma.yy()*sqr(v)
//           + 12.0*sigma.yz()*sigma.yy()*v
//           + w*pow4(v)
//           + 4.0*sigma.yz()*pow3(v)
//         );
// }
//
// IFuncHeader(0,5,0)
// {
//     Is(0,5)[celli] =
//         m0*(15.0*v*sqr(sigma.yy()) + 10.0*sigma.yy()*pow3(v) + pow5(v));
// }
//
// IFuncHeader(1,0,4)
// {
//     Is(1,0,4)[celli] =
//         m0
//        *(
//             3.0*u*sqr(sigma.zz())
//           + 6.0*u*sigma.zz()*sqr(w)
//           + 12.0*sigma.xz()*sigma.zz()*w
//           + u*pow4(w)
//           + 4.0*sigma.xz()*pow3(w)
//         );
// }
//
// IFuncHeader(1,4,0)
// {
//     Is(1,4)[celli] =
//         m0
//        *(
//             3.0*u*sqr(sigma.yy())
//           + 6.0*u*sigma.yy()*sqr(v)
//           + 12.0*sigma.xy()*sigma.yy()*v
//           + u*pow4(v)
//           + 4.0*sigma.xy()*pow3(v)
//         );
// }
//
// IFuncHeader(2,0,3)
// {
//     Is(2,0,3)[celli] =
//         m0
//        *(
//             sqr(u)*pow3(w)
//           + sigma.xx()*pow3(w)
//           + 6.0*sqr(sigma.xz())*w
//           + 6.0*sigma.xz()*u*sqr(w)
//           + 3.0*sigma.zz()*sqr(u)*w
//           + 6.0*sigma.xz()*sigma.zz()*u
//           + 3.0*sigma.xx()*sigma.zz()*w
//         );
// }
//
// IFuncHeader(2,3,0)
// {
//     Is(2,3)[celli] =
//         m0
//        *(
//             sqr(u)*pow3(v)
//           + sigma.xx()*pow3(v)
//           + 6.0*sqr(sigma.xy())*v
//           + 6.0*sigma.xy()*u*sqr(v)
//           + 3.0*sigma.yy()*sqr(u)*v
//           + 6.0*sigma.xy()*sigma.yy()*u
//           + 3.0*sigma.xx()*sigma.yy()*v
//         );
// }
//
// IFuncHeader(3,0,2)
// {
//     Is(3,0,2)[celli] =
//         m0
//        *(
//             pow3(u)*sqr(w)
//           + 6.0*sqr(sigma.xz())*u
//           + sigma.zz()*pow3(u)
//           + 3.0*sigma.xx()*u*sqr(w)
//           + 6.0*sigma.xz()*sqr(u)*w
//           + 3.0*sigma.xx()*sigma.zz()*u
//           + 6.0*sigma.xx()*sigma.xz()*w
//         );
// }
//
// IFuncHeader(3,2,0)
// {
//     Is(3,2)[celli] =
//         m0
//        *(
//             pow3(u)*sqr(v)
//           + 6.0*sqr(sigma.xy())*u
//           + sigma.yy()*pow3(u)
//           + 3.0*sigma.xx()*u*sqr(v)
//           + 6.0*sigma.xy()*sqr(u)*v
//           + 3.0*sigma.xx()*sigma.yy()*u
//           + 6.0*sigma.xx()*sigma.xy()*v
//         );
// }
//
// IFuncHeader(4,0,1)
// {
//     Is(4,0,1)[celli] =
//         m0
//        *(
//             3.0*w*sqr(sigma.xx())
//           + 6.0*w*sigma.xx()*sqr(u)
//           + 12.0*sigma.xz()*sigma.xx()*u
//           + w*pow4(u)
//           + 4.0*sigma.xz()*pow3(u)
//         );
// }
//
// IFuncHeader(4,1,0)
// {
//     Is(4,1)[celli] =
//         m0
//        *(
//             3.0*v*sqr(sigma.xx())
//           + 6.0*v*sigma.xx()*sqr(u)
//           + 12.0*sigma.xy()*sigma.xx()*u
//           + v*pow4(u)
//           + 4.0*sigma.xy()*pow3(u)
//         );
// }
//
// IFuncHeader(5,0,0)
// {
//     Is(5)[celli] =
//         m0*(15.0*u*sqr(sigma.xx()) + 10.0*sigma.xx()*pow3(u) + pow5(u));
// }
//
//
// // Sixth order
// IFuncHeader(0,1,5)
// {
//     Is(0,1,5)[celli] =
//         m0
//        *(
//             15.0*v*w*sqr(sigma.zz())
//           + 15.0*sigma.yz()*sqr(sigma.zz())
//           + 10.0*v*sigma.zz()*pow3(w)
//           + 30.0*sigma.yz()*sigma.zz()*sqr(w)
//           + v*pow5(w)
//           + 5.0*sigma.yz()*pow4(w)
//         );
// }
//
// IFuncHeader(0,2,4)
// {
//     Is(0,2,4)[celli] =
//         m0
//        *(
//             12.0*sqr(sigma.yz())*sigma.zz()
//           + 12.0*sqr(sigma.yz()*w)
//           + 24.0*sigma.yz()*sigma.zz()*v*w
//           + 8.0*sigma.yz()*v*pow3(w)
//           + 3.0*sqr(sigma.zz()*v)
//           + 3.0*sigma.yy()*sqr(sigma.zz())
//           + 6.0*sigma.zz()*sqr(v*w)
//           + 6.0*sigma.yy()*sigma.zz()*sqr(w)
//           + sqr(v*sqr(w))
//           + sigma.yy()*pow4(w)
//         );
// }
//
// IFuncHeader(0,4,2)
// {
//     Is(0,4,2)[celli] =
//         m0
//        *(
//             3.0*sqr(sigma.yy()*w)
//           + 3.0*sigma.zz()*sqr(sigma.yy())
//           + 12.0*sigma.yy()*sqr(sigma.yz())
//           + 24.0*sigma.yy()*sigma.yz()*v*w
//           + 6.0*sigma.yy()*sqr(v*w)
//           + 6.0*sigma.zz()*sigma.yy()*sqr(v)
//           + 12.0*sqr(sigma.yz()*v)
//           + 8.0*sigma.yz()*pow3(v)*w
//           + sqr(sqr(v)*w)
//           + sigma.zz()*pow4(v)
//         );
// }
//
// IFuncHeader(0,5,1)
// {
//     Is(0,5,1)[celli] =
//         m0
//        *(
//             15.0*w*v*sqr(sigma.yy())
//           + 15.0*sigma.yz()*sqr(sigma.yy())
//           + 10.0*w*sigma.yy()*pow3(v)
//           + 30.0*sigma.yz()*sigma.yy()*sqr(v)
//           + w*pow5(v)
//           + 5.0*sigma.yz()*pow4(v)
//         );
// }
//
// IFuncHeader(1,0,5)
// {
//     Is(1,0,5)[celli] =
//         m0
//        *(
//             15.0*u*w*sqr(sigma.zz())
//           + 15.0*sigma.xz()*sqr(sigma.zz())
//           + 10.0*u*sigma.zz()*pow3(w)
//           + 30.0*sigma.xz()*sigma.zz()*sqr(w)
//           + u*pow5(w)
//           + 5.0*sigma.xz()*pow4(w)
//         );
// }
//
// IFuncHeader(1,5,0)
// {
//     Is(1,5)[celli] =
//         m0
//        *(
//             15.0*u*v*sqr(sigma.yy())
//           + 15.0*sigma.xy()*sqr(sigma.yy())
//           + 10.0*u*sigma.yy()*pow3(v)
//           + 30.0*sigma.xy()*sigma.yy()*sqr(v)
//           + u*pow5(v)
//           + 5.0*sigma.xy()*pow4(v)
//         );
// }
//
// IFuncHeader(2,0,4)
// {
//     Is(2,0,4)[celli] =
//         m0
//        *(
//             3.0*sqr(sigma.zz()*u)
//           + 12.0*sqr(sigma.xz()*w)
//           + sqr(u)*pow4(w)
//           + 3.0*sigma.xx()*sqr(sigma.zz())
//           + 12.0*sqr(sigma.xz())*sigma.zz()
//           + sigma.xx()*pow4(w)
//           + 6.0*sigma.xx()*sigma.zz()*sqr(w)
//           + 8.0*sigma.xz()*u*pow3(w)
//           + 6.0*sigma.zz()*sqr(u*w)
//           + 24.0*sigma.xz()*sigma.zz()*u*w
//         );
// }
//
// IFuncHeader(2,4,0)
// {
//     Is(2,4)[celli] =
//         m0
//        *(
//             3.0*sqr(sigma.yy()*u)
//           + 12.0*sqr(sigma.xy()*v)
//           + sqr(u)*pow4(v)
//           + 3.0*sigma.xx()*sqr(sigma.yy())
//           + 12.0*sqr(sigma.xy())*sigma.yy()
//           + sigma.xx()*pow4(v)
//           + 6.0*sigma.xx()*sigma.yy()*sqr(v)
//           + 8.0*sigma.xy()*u*pow3(v)
//           + 6.0*sigma.yy()*sqr(u*v)
//           + 24.0*sigma.xy()*sigma.yy()*u*v
//         );
// }
//
// IFuncHeader(4,0,2)
// {
//     Is(4,0,2)[celli] =
//         m0
//        *(
//             12.0*sqr(sigma.xz()*u)
//           + 3.0*sqr(sigma.xx()*w)
//           + pow4(u)*sqr(w)
//           + 12.0*sigma.xx()*sqr(sigma.xz())
//           + 3.0*sqr(sigma.xx())*sigma.zz()
//           + sigma.zz()*pow4(u)
//           + 6.0*sigma.xx()*sigma.zz()*sqr(u)
//           + 8.0*sigma.xz()*pow3(u)*w
//           + 6.0*sigma.xx()*sqr(u*w)
//           + 24.0*sigma.xx()*sigma.xz()*u*w
//         );
// }
//
// IFuncHeader(4,2,0)
// {
//     Is(4,2)[celli] =
//         m0
//        *(
//             12.0*sqr(sigma.xy()*u)
//           + 3.0*sqr(sigma.xx()*v)
//           + pow4(u)*sqr(v)
//           + 12.0*sigma.xx()*sqr(sigma.xy())
//           + 3.0*sqr(sigma.xx())*sigma.yy()
//           + sigma.yy()*pow4(u)
//           + 6.0*sigma.xx()*sigma.yy()*sqr(u)
//           + 8.0*sigma.xy()*pow3(u)*v
//           + 6.0*sigma.xx()*sqr(u*v)
//           + 24.0*sigma.xx()*sigma.xy()*u*v
//         );
// }
//
// IFuncHeader(5,0,1)
// {
//     Is(5,0,1)[celli] =
//         m0
//        *(
//             15.0*sqr(sigma.xx())*sigma.xz()
//           + 5.0*sigma.xz()*pow4(u)
//           + pow5(u)*w
//           + 30.0*sigma.xx()*sigma.xz()*sqr(u)
//           + 15.0*sqr(sigma.xx())*u*w
//           + 10.0*sigma.xx()*pow3(u)*w
//         );
// }
//
// IFuncHeader(5,1,0)
// {
//     Is(5,1)[celli] =
//         m0
//        *(
//             15.0*sqr(sigma.xx())*sigma.xy()
//           + 5.0*sigma.xy()*pow4(u)
//           + pow5(u)*v
//           + 30.0*sigma.xx()*sigma.xy()*sqr(u)
//           + 15.0*sqr(sigma.xx())*u*v
//           + 10.0*sigma.xx()*pow3(u)*v
//         );
// }
//
// IFuncHeader(0,2,5)
// {
//     Is(0,2,5)[celli] =
//         m0
//        *(
//             20.0*sqr(sigma.yz())*pow3(w)
//           + sqr(v)*pow5(w)
//           + sigma.yy()*pow5(w)
//           + 30.0*sigma.yz()*sqr(sigma.zz())*v
//           + 15.0*sigma.yy()*sqr(sigma.zz())*w
//           + 10.0*sigma.yy()*sigma.zz()*pow3(w)
//           + 60.0*sqr(sigma.yz())*sigma.zz()*w
//           + 10.0*sigma.yz()*v*pow4(w)
//           + 15.0*sqr(sigma.zz())*sqr(v)*w
//           + 10.0*sigma.zz()*sqr(v)*pow3(w)
//           + 60.0*sigma.yz()*sigma.zz()*v*sqr(w)
//         );
// }
//
// IFuncHeader(0,5,2)
// {
//     Is(0,5,2)[celli] =
//         m0
//        *(
//             20.0*sqr(sigma.yz())*pow3(v)
//           + pow5(v)*sqr(w)
//           + sigma.zz()*pow5(v)
//           + 60.0*sigma.yy()*sqr(sigma.yz())*v
//           + 15.0*sqr(sigma.yy())*sigma.zz()*v
//           + 10.0*sigma.yy()*sigma.zz()*pow3(v)
//           + 30.0*sqr(sigma.yy())*sigma.yz()*w
//           + 10.0*sigma.yz()*pow4(v)*w
//           + 15.0*sqr(sigma.yy())*v*sqr(w)
//           + 10.0*sigma.yy()*pow3(v)*sqr(w)
//           + 60.0*sigma.yy()*sigma.yz()*sqr(v)*w
//         );
// }
//
// IFuncHeader(2,0,5)
// {
//     Is(2,0,5)[celli] =
//         m0
//        *(
//             20.0*sqr(sigma.xz())*pow3(w)
//           + sqr(u)*pow5(w)
//           + sigma.xx()*pow5(w)
//           + 30.0*sigma.xz()*sqr(sigma.zz())*u
//           + 15.0*sigma.xx()*sqr(sigma.zz())*w
//           + 10.0*sigma.xx()*sigma.zz()*pow3(w)
//           + 60.0*sqr(sigma.xz())*sigma.zz()*w
//           + 10.0*sigma.xz()*u*pow4(w)
//           + 15.0*sqr(sigma.zz())*sqr(u)*w
//           + 10.0*sigma.zz()*sqr(u)*pow3(w)
//           + 60.0*sigma.xz()*sigma.zz()*u*sqr(w)
//         );
// }
//
// IFuncHeader(2,5,0)
// {
//     Is(2,5)[celli] =
//         m0
//        *(
//             20.0*sqr(sigma.xy())*pow3(v)
//           + sqr(u)*pow5(v)
//           + sigma.xx()*pow5(v)
//           + 30.0*sigma.xy()*sqr(sigma.yy())*u
//           + 15.0*sigma.xx()*sqr(sigma.yy())*v
//           + 10.0*sigma.xx()*sigma.yy()*pow3(v)
//           + 60.0*sqr(sigma.xy())*sigma.yy()*v
//           + 10.0*sigma.xy()*u*pow4(v)
//           + 15.0*sqr(sigma.yy())*sqr(u)*v
//           + 10.0*sigma.yy()*sqr(u)*pow3(v)
//           + 60.0*sigma.xy()*sigma.yy()*u*sqr(v)
//         );
// }
//
// IFuncHeader(5,0,2)
// {
//     Is(5,0,2)[celli] =
//         m0
//        *(
//             20.0*sqr(sigma.xz())*pow3(u)
//           + pow5(u)*sqr(w)
//           + sigma.zz()*pow5(u)
//           + 60.0*sigma.xx()*sqr(sigma.xz())*u
//           + 15.0*sqr(sigma.xx())*sigma.zz()*u
//           + 10.0*sigma.xx()*sigma.zz()*pow3(u)
//           + 30.0*sqr(sigma.xx())*sigma.xz()*w
//           + 10.0*sigma.xz()*pow4(u)*w
//           + 15.0*sqr(sigma.xx())*u*sqr(w)
//           + 10.0*sigma.xx()*pow3(u)*sqr(w)
//           + 60.0*sigma.xx()*sigma.xz()*sqr(u)*w
//         );
// }
//
// IFuncHeader(5,2,0)
// {
//     Is(5,2)[celli] =
//         m0
//        *(
//             20.0*sqr(sigma.xy())*pow3(u)
//           + pow5(u)*sqr(v)
//           + sigma.yy()*pow5(u)
//           + 60.0*sigma.xx()*sqr(sigma.xy())*u
//           + 15.0*sqr(sigma.xx())*sigma.yy()*u
//           + 10.0*sigma.xx()*sigma.yy()*pow3(u)
//           + 30.0*sqr(sigma.xx())*sigma.xy()*v
//           + 10.0*sigma.xy()*pow4(u)*v
//           + 15.0*sqr(sigma.xx())*u*sqr(v)
//           + 10.0*sigma.xx()*pow3(u)*sqr(v)
//           + 60.0*sigma.xx()*sigma.xy()*sqr(u)*v
//         );
// }
// ************************************************************************* //
