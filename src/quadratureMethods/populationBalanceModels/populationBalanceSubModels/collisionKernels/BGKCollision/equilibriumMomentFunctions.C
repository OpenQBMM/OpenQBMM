/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | OpenQBMM - www.openqbmm.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Code created 2015-2018 by Alberto Passalacqua
    Contributed 2018-07-31 to the OpenFOAM Foundation
    Copyright (C) 2018 OpenFOAM Foundation
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
#include "BGKCollision.H"

// Zero order
momentFuncHeader(0,0,0)
{
    moments(0) = m0;
}

// First order
momentFuncHeader(0,0,1)
{
    moments(0,0,1) = m0*w;
}

momentFuncHeader(0,1,0)
{
    moments(0,1) = m0*v;
}

momentFuncHeader(1,0,0)
{
    moments(1) = m0*u;
}

// Second order
momentFuncHeader(0,0,2)
{
    moments(0,0,2) = m0*(sigma.zz() + sqr(w));
}

momentFuncHeader(0,1,1)
{
    moments(0,1,1) = m0*(sigma.yz() + v*w);
}

momentFuncHeader(1,0,1)
{
    moments(1,0,1) = m0*(sigma.xz() + u*w);
}

momentFuncHeader(1,1,0)
{
    moments(1,1) = m0*(sigma.xy() + u*v);
}

momentFuncHeader(0,2,0)
{
    moments(0,2) = m0*(sigma.yy() + sqr(v));
}

momentFuncHeader(2,0,0)
{
    moments(2) = m0*(sigma.xx() + sqr(u));
}

// Third order
momentFuncHeader(0,0,3)
{
    moments(0,0,3) = m0*(3.0*sigma.zz()*w + pow3(w));
}

momentFuncHeader(0,1,2)
{
    moments(0,1,2) = m0*(v*sqr(w) + 2.0*sigma.yz()*w + sigma.zz()*v);
}

momentFuncHeader(0,2,1)
{
    moments(0,2,1) = m0*(w*sqr(v) + 2.0*sigma.yz()*v + sigma.yy()*w);
}

momentFuncHeader(0,3,0)
{
    moments(0,3) = m0*(3.0*sigma.yy()*v + pow3(v));
}

momentFuncHeader(1,0,2)
{
    moments(1,0,2) = m0*(u*sqr(w) + 2.0*sigma.xz()*w + sigma.zz()*u);
}

momentFuncHeader(1,1,1)
{
    moments(1,1,1) = m0*(sigma.xy()*w + sigma.xz()*v + sigma.yz()*u + u*v*w);
}

momentFuncHeader(1,2,0)
{
    moments(1,2) = m0*(u*sqr(v) + 2.0*sigma.xy()*v + sigma.yy()*u);
}

momentFuncHeader(2,0,1)
{
    moments(2,0,1) = m0*(w*sqr(u) + 2.0*sigma.xz()*u + sigma.xx()*w);
}

momentFuncHeader(2,1,0)
{
    moments(2,1) = m0*(v*sqr(u) + 2.0*sigma.xy()*u + sigma.xx()*v);
}

momentFuncHeader(3,0,0)
{
    moments(3) = m0*(3.0*sigma.xx()*u + pow3(u));
}

// Fourth order
momentFuncHeader(0,0,4)
{
    moments(0,0,4) =
        m0*(6.0*sqr(w)*sigma.zz() + 3.0*sqr(sigma.zz()) + pow4(w));
}

momentFuncHeader(0,1,3)
{
    moments(0,1,3) =
        m0
       *(
            v*pow3(w)
          + 3.0*sigma.yz()*sqr(w)
          + 3.0*sigma.zz()*w*v
          + 3.0*sigma.zz()*sigma.yz()
        );
}

momentFuncHeader(0,2,2)
{
    moments(0,2,2) =
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

momentFuncHeader(0,3,1)
{
    moments(0,3,1) =
        m0
       *(
            w*pow3(v)
          + 3.0*sigma.yz()*sqr(v)
          + 3.0*sigma.yy()*v*w
          + 3.0*sigma.yy()*sigma.yz()
        );
}

momentFuncHeader(0,4,0)
{
    moments(0,4) =
        m0*(6.0*sqr(v)*sigma.yy() + 3.0*sqr(sigma.yy()) + pow4(v));
}

momentFuncHeader(1,0,3)
{
    moments(1,0,3) =
        m0
       *(
            u*pow3(w)
          + 3.0*sigma.xz()*sqr(w)
          + 3.0*sigma.zz()*w*u
          + 3.0*sigma.zz()*sigma.xz()
        );
}

momentFuncHeader(2,0,2)
{
    moments(2,0,2) =
        m0
       *(
            2.0*sqr(sigma.xz())
          + 4.0*sigma.xz()*u*w
          + sqr(u*w)
          + sigma.zz()*sqr(u) + sigma.xx()*sqr(w)
          + sigma.xx()*sigma.zz()
        );
}

momentFuncHeader(2,2,0)
{
    moments(2,2) =
        m0
       *(
            2.0*sqr(sigma.xy())
          + 4.0*sigma.xy()*u*v
          + sqr(u*v)
          + sigma.yy()*sqr(u)
          + sigma.xx()*sqr(v)
          + sigma.xx()*sigma.yy()
        );
}

momentFuncHeader(1,3,0)
{
    moments(1,3) =
        m0
       *(
            u*pow3(v)
          + 3.0*sigma.xy()*sqr(v)
          + 3.0*sigma.yy()*v*u
          + 3.0*sigma.yy()*sigma.xy()
        );
}

momentFuncHeader(3,0,1)
{
    moments(3,0,1) =
        m0
       *(
            w*pow3(u)
          + 3.0*sigma.xz()*sqr(u)
          + 3.0*sigma.xx()*u*w
          + 3.0*sigma.xx()*sigma.xz()
        );
}

momentFuncHeader(3,1,0)
{
    moments(3,1) =
        m0
       *(
            v*pow3(u)
          + 3.0*sigma.xy()*sqr(u)
          + 3.0*sigma.xx()*u*v
          + 3.0*sigma.xx()*sigma.xy()
        );
}

momentFuncHeader(4,0,0)
{
    moments(4) =
        m0*(6.0*sqr(u)*sigma.xx() + 3.0*sqr(sigma.xx()) + pow4(u));
}


// Fith order
momentFuncHeader(0,0,5)
{
    moments(0,0,5) =
        m0*(15.0*w*sqr(sigma.zz()) + 10.0*sigma.zz()*pow3(w) + pow5(w));
}

momentFuncHeader(0,1,4)
{
    moments(0,1,4) =
        m0
       *(
            3.0*v*sqr(sigma.zz())
          + 6.0*v*sigma.zz()*sqr(w)
          + 12.0*sigma.yz()*sigma.zz()*w
          + v*pow4(w)
          + 4.0*sigma.yz()*pow3(w)
        );
}

momentFuncHeader(0,2,3)
{
    moments(0,2,3) =
        m0
       *(
            sqr(v)*pow3(w)
          + sigma.yy()*pow3(w)
          + 6.0*sqr(sigma.yz())*w
          + 6.0*sigma.yz()*v*sqr(w)
          + 3.0*sigma.zz()*sqr(v)*w
          + 6.0*sigma.yz()*sigma.zz()*v
          + 3.0*sigma.yy()*sigma.zz()*w
        );
}

momentFuncHeader(0,3,2)
{
    moments(0,3,2) =
        m0
       *(
            pow3(v)*sqr(w)
          + 6.0*sqr(sigma.yz())*v
          + sigma.zz()*pow3(v)
          + 3.0*sigma.yy()*v*sqr(w)
          + 6.0*sigma.yz()*sqr(v)*w
          + 3.0*sigma.yy()*sigma.zz()*v
          + 6.0*sigma.yy()*sigma.yz()*w
        );
}

momentFuncHeader(0,4,1)
{
    moments(0,4,1) =
        m0
       *(
            3.0*w*sqr(sigma.yy())
          + 6.0*w*sigma.yy()*sqr(v)
          + 12.0*sigma.yz()*sigma.yy()*v
          + w*pow4(v)
          + 4.0*sigma.yz()*pow3(v)
        );
}

momentFuncHeader(0,5,0)
{
    moments(0,5) =
        m0*(15.0*v*sqr(sigma.yy()) + 10.0*sigma.yy()*pow3(v) + pow5(v));
}

momentFuncHeader(1,0,4)
{
    moments(1,0,4) =
        m0
       *(
            3.0*u*sqr(sigma.zz())
          + 6.0*u*sigma.zz()*sqr(w)
          + 12.0*sigma.xz()*sigma.zz()*w
          + u*pow4(w)
          + 4.0*sigma.xz()*pow3(w)
        );
}

momentFuncHeader(1,4,0)
{
    moments(1,4) =
        m0
       *(
            3.0*u*sqr(sigma.yy())
          + 6.0*u*sigma.yy()*sqr(v)
          + 12.0*sigma.xy()*sigma.yy()*v
          + u*pow4(v)
          + 4.0*sigma.xy()*pow3(v)
        );
}

momentFuncHeader(2,0,3)
{
    moments(2,0,3) =
        m0
       *(
            sqr(u)*pow3(w)
          + sigma.xx()*pow3(w)
          + 6.0*sqr(sigma.xz())*w
          + 6.0*sigma.xz()*u*sqr(w)
          + 3.0*sigma.zz()*sqr(u)*w
          + 6.0*sigma.xz()*sigma.zz()*u
          + 3.0*sigma.xx()*sigma.zz()*w
        );
}

momentFuncHeader(2,3,0)
{
    moments(2,3) =
        m0
       *(
            sqr(u)*pow3(v)
          + sigma.xx()*pow3(v)
          + 6.0*sqr(sigma.xy())*v
          + 6.0*sigma.xy()*u*sqr(v)
          + 3.0*sigma.yy()*sqr(u)*v
          + 6.0*sigma.xy()*sigma.yy()*u
          + 3.0*sigma.xx()*sigma.yy()*v
        );
}

momentFuncHeader(3,0,2)
{
    moments(3,0,2) =
        m0
       *(
            pow3(u)*sqr(w)
          + 6.0*sqr(sigma.xz())*u
          + sigma.zz()*pow3(u)
          + 3.0*sigma.xx()*u*sqr(w)
          + 6.0*sigma.xz()*sqr(u)*w
          + 3.0*sigma.xx()*sigma.zz()*u
          + 6.0*sigma.xx()*sigma.xz()*w
        );
}

momentFuncHeader(3,2,0)
{
    moments(3,2) =
        m0
       *(
            pow3(u)*sqr(v)
          + 6.0*sqr(sigma.xy())*u
          + sigma.yy()*pow3(u)
          + 3.0*sigma.xx()*u*sqr(v)
          + 6.0*sigma.xy()*sqr(u)*v
          + 3.0*sigma.xx()*sigma.yy()*u
          + 6.0*sigma.xx()*sigma.xy()*v
        );
}

momentFuncHeader(4,0,1)
{
    moments(4,0,1) =
        m0
       *(
            3.0*w*sqr(sigma.xx())
          + 6.0*w*sigma.xx()*sqr(u)
          + 12.0*sigma.xz()*sigma.xx()*u
          + w*pow4(u)
          + 4.0*sigma.xz()*pow3(u)
        );
}

momentFuncHeader(4,1,0)
{
    moments(4,1) =
        m0
       *(
            3.0*v*sqr(sigma.xx())
          + 6.0*v*sigma.xx()*sqr(u)
          + 12.0*sigma.xy()*sigma.xx()*u
          + v*pow4(u)
          + 4.0*sigma.xy()*pow3(u)
        );
}

momentFuncHeader(5,0,0)
{
    moments(5) =
        m0*(15.0*u*sqr(sigma.xx()) + 10.0*sigma.xx()*pow3(u) + pow5(u));
}


// Sixth order
momentFuncHeader(0,1,5)
{
    moments(0,1,5) =
        m0
       *(
            15.0*v*w*sqr(sigma.zz())
          + 15.0*sigma.yz()*sqr(sigma.zz())
          + 10.0*v*sigma.zz()*pow3(w)
          + 30.0*sigma.yz()*sigma.zz()*sqr(w)
          + v*pow5(w)
          + 5.0*sigma.yz()*pow4(w)
        );
}

momentFuncHeader(0,2,4)
{
    moments(0,2,4) =
        m0
       *(
            12.0*sqr(sigma.yz())*sigma.zz()
          + 12.0*sqr(sigma.yz()*w)
          + 24.0*sigma.yz()*sigma.zz()*v*w
          + 8.0*sigma.yz()*v*pow3(w)
          + 3.0*sqr(sigma.zz()*v)
          + 3.0*sigma.yy()*sqr(sigma.zz())
          + 6.0*sigma.zz()*sqr(v*w)
          + 6.0*sigma.yy()*sigma.zz()*sqr(w)
          + sqr(v*sqr(w))
          + sigma.yy()*pow4(w)
        );
}

momentFuncHeader(0,4,2)
{
    moments(0,4,2) =
        m0
       *(
            3.0*sqr(sigma.yy()*w)
          + 3.0*sigma.zz()*sqr(sigma.yy())
          + 12.0*sigma.yy()*sqr(sigma.yz())
          + 24.0*sigma.yy()*sigma.yz()*v*w
          + 6.0*sigma.yy()*sqr(v*w)
          + 6.0*sigma.zz()*sigma.yy()*sqr(v)
          + 12.0*sqr(sigma.yz()*v)
          + 8.0*sigma.yz()*pow3(v)*w
          + sqr(sqr(v)*w)
          + sigma.zz()*pow4(v)
        );
}

momentFuncHeader(0,5,1)
{
    moments(0,5,1) =
        m0
       *(
            15.0*w*v*sqr(sigma.yy())
          + 15.0*sigma.yz()*sqr(sigma.yy())
          + 10.0*w*sigma.yy()*pow3(v)
          + 30.0*sigma.yz()*sigma.yy()*sqr(v)
          + w*pow5(v)
          + 5.0*sigma.yz()*pow4(v)
        );
}

momentFuncHeader(1,0,5)
{
    moments(1,0,5) =
        m0
       *(
            15.0*u*w*sqr(sigma.zz())
          + 15.0*sigma.xz()*sqr(sigma.zz())
          + 10.0*u*sigma.zz()*pow3(w)
          + 30.0*sigma.xz()*sigma.zz()*sqr(w)
          + u*pow5(w)
          + 5.0*sigma.xz()*pow4(w)
        );
}

momentFuncHeader(1,5,0)
{
    moments(1,5) =
        m0
       *(
            15.0*u*v*sqr(sigma.yy())
          + 15.0*sigma.xy()*sqr(sigma.yy())
          + 10.0*u*sigma.yy()*pow3(v)
          + 30.0*sigma.xy()*sigma.yy()*sqr(v)
          + u*pow5(v)
          + 5.0*sigma.xy()*pow4(v)
        );
}

momentFuncHeader(2,0,4)
{
    moments(2,0,4) =
        m0
       *(
            3.0*sqr(sigma.zz()*u)
          + 12.0*sqr(sigma.xz()*w)
          + sqr(u)*pow4(w)
          + 3.0*sigma.xx()*sqr(sigma.zz())
          + 12.0*sqr(sigma.xz())*sigma.zz()
          + sigma.xx()*pow4(w)
          + 6.0*sigma.xx()*sigma.zz()*sqr(w)
          + 8.0*sigma.xz()*u*pow3(w)
          + 6.0*sigma.zz()*sqr(u*w)
          + 24.0*sigma.xz()*sigma.zz()*u*w
        );
}

momentFuncHeader(2,4,0)
{
    moments(2,4) =
        m0
       *(
            3.0*sqr(sigma.yy()*u)
          + 12.0*sqr(sigma.xy()*v)
          + sqr(u)*pow4(v)
          + 3.0*sigma.xx()*sqr(sigma.yy())
          + 12.0*sqr(sigma.xy())*sigma.yy()
          + sigma.xx()*pow4(v)
          + 6.0*sigma.xx()*sigma.yy()*sqr(v)
          + 8.0*sigma.xy()*u*pow3(v)
          + 6.0*sigma.yy()*sqr(u*v)
          + 24.0*sigma.xy()*sigma.yy()*u*v
        );
}

momentFuncHeader(4,0,2)
{
    moments(4,0,2) =
        m0
       *(
            12.0*sqr(sigma.xz()*u)
          + 3.0*sqr(sigma.xx()*w)
          + pow4(u)*sqr(w)
          + 12.0*sigma.xx()*sqr(sigma.xz())
          + 3.0*sqr(sigma.xx())*sigma.zz()
          + sigma.zz()*pow4(u)
          + 6.0*sigma.xx()*sigma.zz()*sqr(u)
          + 8.0*sigma.xz()*pow3(u)*w
          + 6.0*sigma.xx()*sqr(u*w)
          + 24.0*sigma.xx()*sigma.xz()*u*w
        );
}

momentFuncHeader(4,2,0)
{
    moments(4,2) =
        m0
       *(
            12.0*sqr(sigma.xy()*u)
          + 3.0*sqr(sigma.xx()*v)
          + pow4(u)*sqr(v)
          + 12.0*sigma.xx()*sqr(sigma.xy())
          + 3.0*sqr(sigma.xx())*sigma.yy()
          + sigma.yy()*pow4(u)
          + 6.0*sigma.xx()*sigma.yy()*sqr(u)
          + 8.0*sigma.xy()*pow3(u)*v
          + 6.0*sigma.xx()*sqr(u*v)
          + 24.0*sigma.xx()*sigma.xy()*u*v
        );
}

momentFuncHeader(5,0,1)
{
    moments(5,0,1) =
        m0
       *(
            15.0*sqr(sigma.xx())*sigma.xz()
          + 5.0*sigma.xz()*pow4(u)
          + pow5(u)*w
          + 30.0*sigma.xx()*sigma.xz()*sqr(u)
          + 15.0*sqr(sigma.xx())*u*w
          + 10.0*sigma.xx()*pow3(u)*w
        );
}

momentFuncHeader(5,1,0)
{
    moments(5,1) =
        m0
       *(
            15.0*sqr(sigma.xx())*sigma.xy()
          + 5.0*sigma.xy()*pow4(u)
          + pow5(u)*v
          + 30.0*sigma.xx()*sigma.xy()*sqr(u)
          + 15.0*sqr(sigma.xx())*u*v
          + 10.0*sigma.xx()*pow3(u)*v
        );
}

momentFuncHeader(0,2,5)
{
    moments(0,2,5) =
        m0
       *(
            20.0*sqr(sigma.yz())*pow3(w)
          + sqr(v)*pow5(w)
          + sigma.yy()*pow5(w)
          + 30.0*sigma.yz()*sqr(sigma.zz())*v
          + 15.0*sigma.yy()*sqr(sigma.zz())*w
          + 10.0*sigma.yy()*sigma.zz()*pow3(w)
          + 60.0*sqr(sigma.yz())*sigma.zz()*w
          + 10.0*sigma.yz()*v*pow4(w)
          + 15.0*sqr(sigma.zz())*sqr(v)*w
          + 10.0*sigma.zz()*sqr(v)*pow3(w)
          + 60.0*sigma.yz()*sigma.zz()*v*sqr(w)
        );
}

momentFuncHeader(0,5,2)
{
    moments(0,5,2) =
        m0
       *(
            20.0*sqr(sigma.yz())*pow3(v)
          + pow5(v)*sqr(w)
          + sigma.zz()*pow5(v)
          + 60.0*sigma.yy()*sqr(sigma.yz())*v
          + 15.0*sqr(sigma.yy())*sigma.zz()*v
          + 10.0*sigma.yy()*sigma.zz()*pow3(v)
          + 30.0*sqr(sigma.yy())*sigma.yz()*w
          + 10.0*sigma.yz()*pow4(v)*w
          + 15.0*sqr(sigma.yy())*v*sqr(w)
          + 10.0*sigma.yy()*pow3(v)*sqr(w)
          + 60.0*sigma.yy()*sigma.yz()*sqr(v)*w
        );
}

momentFuncHeader(2,0,5)
{
    moments(2,0,5) =
        m0
       *(
            20.0*sqr(sigma.xz())*pow3(w)
          + sqr(u)*pow5(w)
          + sigma.xx()*pow5(w)
          + 30.0*sigma.xz()*sqr(sigma.zz())*u
          + 15.0*sigma.xx()*sqr(sigma.zz())*w
          + 10.0*sigma.xx()*sigma.zz()*pow3(w)
          + 60.0*sqr(sigma.xz())*sigma.zz()*w
          + 10.0*sigma.xz()*u*pow4(w)
          + 15.0*sqr(sigma.zz())*sqr(u)*w
          + 10.0*sigma.zz()*sqr(u)*pow3(w)
          + 60.0*sigma.xz()*sigma.zz()*u*sqr(w)
        );
}

momentFuncHeader(2,5,0)
{
    moments(2,5) =
        m0
       *(
            20.0*sqr(sigma.xy())*pow3(v)
          + sqr(u)*pow5(v)
          + sigma.xx()*pow5(v)
          + 30.0*sigma.xy()*sqr(sigma.yy())*u
          + 15.0*sigma.xx()*sqr(sigma.yy())*v
          + 10.0*sigma.xx()*sigma.yy()*pow3(v)
          + 60.0*sqr(sigma.xy())*sigma.yy()*v
          + 10.0*sigma.xy()*u*pow4(v)
          + 15.0*sqr(sigma.yy())*sqr(u)*v
          + 10.0*sigma.yy()*sqr(u)*pow3(v)
          + 60.0*sigma.xy()*sigma.yy()*u*sqr(v)
        );
}

momentFuncHeader(5,0,2)
{
    moments(5,0,2) =
        m0
       *(
            20.0*sqr(sigma.xz())*pow3(u)
          + pow5(u)*sqr(w)
          + sigma.zz()*pow5(u)
          + 60.0*sigma.xx()*sqr(sigma.xz())*u
          + 15.0*sqr(sigma.xx())*sigma.zz()*u
          + 10.0*sigma.xx()*sigma.zz()*pow3(u)
          + 30.0*sqr(sigma.xx())*sigma.xz()*w
          + 10.0*sigma.xz()*pow4(u)*w
          + 15.0*sqr(sigma.xx())*u*sqr(w)
          + 10.0*sigma.xx()*pow3(u)*sqr(w)
          + 60.0*sigma.xx()*sigma.xz()*sqr(u)*w
        );
}

momentFuncHeader(5,2,0)
{
    moments(5,2) =
        m0
       *(
            20.0*sqr(sigma.xy())*pow3(u)
          + pow5(u)*sqr(v)
          + sigma.yy()*pow5(u)
          + 60.0*sigma.xx()*sqr(sigma.xy())*u
          + 15.0*sqr(sigma.xx())*sigma.yy()*u
          + 10.0*sigma.xx()*sigma.yy()*pow3(u)
          + 30.0*sqr(sigma.xx())*sigma.xy()*v
          + 10.0*sigma.xy()*pow4(u)*v
          + 15.0*sqr(sigma.xx())*u*sqr(v)
          + 10.0*sigma.xx()*pow3(u)*sqr(v)
          + 60.0*sigma.xx()*sigma.xy()*sqr(u)*v
        );
}
// ************************************************************************* //
