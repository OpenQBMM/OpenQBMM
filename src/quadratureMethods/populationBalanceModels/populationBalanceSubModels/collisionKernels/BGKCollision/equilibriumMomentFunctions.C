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

momentFuncHeader(0,0,0)
{
    moments(0)[celli] = m0;
}

momentFuncHeader(1,0,0)
{
    moments(1)[celli] = m0*u;
}

momentFuncHeader(0,1,0)
{
    moments(0,1)[celli] = m0*v;
}

momentFuncHeader(0,0,1)
{
    moments(0,0,1)[celli] = m0*w;
}

momentFuncHeader(2,0,0)
{
    moments(2)[celli] = m0*(sigma.xx() + sqr(u));
}

momentFuncHeader(0,2,0)
{
    moments(0,2)[celli] = m0*(sigma.yy() + sqr(v));
}

momentFuncHeader(0,0,2)
{
    moments(0,0,2)[celli] = m0*(sigma.zz() + sqr(w));
}

momentFuncHeader(3,0,0)
{
    moments(3)[celli] = m0*(3.0*sigma.xx()*u + pow3(u));
}

momentFuncHeader(0,3,0)
{
    moments(0,3)[celli] = m0*(3.0*sigma.yy()*v + pow3(v));
}

momentFuncHeader(0,0,3)
{
    moments(0,0,3)[celli] = m0*(3.0*sigma.zz()*w + pow3(w));
}

momentFuncHeader(4,0,0)
{
    moments(4)[celli] =
        m0*(6.0*sqr(u)*sigma.xx() + 3.0*sqr(sigma.xx()) + pow4(u));
}

momentFuncHeader(0,4,0)
{
    moments(0,4)[celli] =
        m0*(6.0*sqr(v)*sigma.yy() + 3.0*sqr(sigma.yy()) + pow4(v));
}

momentFuncHeader(0,0,4)
{
    moments(0,0,4)[celli] =
        m0*(6.0*sqr(w)*sigma.zz() + 3.0*sqr(sigma.zz()) + pow4(w));
}

momentFuncHeader(5,0,0)
{
    moments(5)[celli] =
        m0*(15.0*u*sqr(sigma.xx()) + 10.0*sigma.xx()*pow3(u) + pow5(u));
}

momentFuncHeader(0,5,0)
{
    moments(0,5)[celli] =
        m0*(15.0*v*sqr(sigma.yy()) + 10.0*sigma.yy()*pow3(v) + pow5(v));
}

momentFuncHeader(0,0,5)
{
    moments(0,0,5)[celli] =
        m0*(15.0*w*sqr(sigma.zz()) + 10.0*sigma.zz()*pow3(w) + pow5(w));
}

momentFuncHeader(1,1,0)
{
    moments(1,1,0)[celli] = m0*(sigma.xy() + u*v);
}

momentFuncHeader(1,0,1)
{
    moments(1,0,1)[celli] = m0*(sigma.xz() + u*w);
}

momentFuncHeader(0,1,1)
{
    moments(0,1,1)[celli] = m0*(sigma.yz() + v*w);
}

momentFuncHeader(2,1,0)
{
    moments(2,1,0)[celli] = m0*(v*sqr(u) + 2.0*sigma.xy()*u + sigma.xx()*v);
}

momentFuncHeader(2,0,1)
{
    moments(2,0,1)[celli] = m0*(w*sqr(u) + 2.0*sigma.xz()*u + sigma.xx()*w);
}

momentFuncHeader(0,2,1)
{
    moments(0,2,1)[celli] = m0*(w*sqr(v) + 2.0*sigma.yz()*v + sigma.yy()*w);
}

momentFuncHeader(1,2,0)
{
    moments(1,2,0)[celli] = m0*(u*sqr(v) + 2.0*sigma.xy()*v + sigma.yy()*u);
}

momentFuncHeader(1,0,2)
{
    moments(1,0,2)[celli] = m0*(u*sqr(w) + 2.0*sigma.xz()*w + sigma.zz()*u);
}

momentFuncHeader(0,1,2)
{
    moments(0,1,2)[celli] = m0*(v*sqr(w) + 2.0*sigma.yz()*w + sigma.zz()*v);
}

momentFuncHeader(3,1,0)
{
    moments(3,1,0)[celli] =
        m0
       *(
            v*pow3(u)
          + 3.0*sigma.xy()*sqr(u)
          + 3.0*sigma.xx()*u*v
          + 3.0*sigma.xx()*sigma.xy()
        );
}

momentFuncHeader(3,0,1)
{
    moments(3,0,1)[celli] =
        m0
       *(
            w*pow3(u)
          + 3.0*sigma.xz()*sqr(u)
          + 3.0*sigma.xx()*u*w
          + 3.0*sigma.xx()*sigma.xz()
        );
}

momentFuncHeader(0,3,1)
{
    moments(0,3,1)[celli] =
        m0
       *(
            w*pow3(v)
          + 3.0*sigma.yz()*sqr(v)
          + 3.0*sigma.yy()*v*w
          + 3.0*sigma.yy()*sigma.yz()
        );
}

momentFuncHeader(1,3,0)
{
    moments(1,3,0)[celli] =
        m0
       *(
            u*pow3(v)
          + 3.0*sigma.xy()*sqr(v)
          + 3.0*sigma.yy()*v*u
          + 3.0*sigma.yy()*sigma.xy()
        );
}

momentFuncHeader(1,0,3)
{
    moments(1,0,3)[celli] =
        m0
       *(
            u*pow3(w)
          + 3.0*sigma.xz()*sqr(w)
          + 3.0*sigma.zz()*w*u
          + 3.0*sigma.zz()*sigma.xz()
        );
}

momentFuncHeader(0,1,3)
{
    moments(0,1,3)[celli] =
        m0
       *(
            v*pow3(w)
          + 3.0*sigma.yz()*sqr(w)
          + 3.0*sigma.zz()*w*v
          + 3.0*sigma.zz()*sigma.yz()
        );
}

momentFuncHeader(4,1,0)
{
    moments(4,1,0)[celli] =
        m0
       *(
            3.0*v*sqr(sigma.xx())
          + 6.0*v*sigma.xx()*sqr(u)
          + 12.0*sigma.xy()*sigma.xx()*u
          + v*pow4(u)
          + 4.0*sigma.xy()*pow3(u)
        );
}

momentFuncHeader(4,0,1)
{
    moments(4,0,1)[celli] =
        m0
       *(
            3.0*w*sqr(sigma.xx())
          + 6.0*w*sigma.xx()*sqr(u)
          + 12.0*sigma.xz()*sigma.xx()*u
          + w*pow4(u)
          + 4.0*sigma.xz()*pow3(u)
        );
}

momentFuncHeader(0,4,1)
{
    moments(0,4,1)[celli] =
        m0
       *(
            3.0*w*sqr(sigma.yy())
          + 6.0*w*sigma.yy()*sqr(v)
          + 12.0*sigma.yz()*sigma.yy()*v
          + w*pow4(v)
          + 4.0*sigma.yz()*pow3(v)
        );
}

momentFuncHeader(1,4,0)
{
    moments(1,4,0)[celli] =
        m0
       *(
            3.0*u*sqr(sigma.yy())
          + 6.0*u*sigma.yy()*sqr(v)
          + 12.0*sigma.xy()*sigma.yy()*v
          + u*pow4(v)
          + 4.0*sigma.xy()*pow3(v)
        );
}

momentFuncHeader(1,0,4)
{
    moments(1,0,4)[celli] =
        m0
       *(
            3.0*u*sqr(sigma.zz())
          + 6.0*u*sigma.zz()*sqr(w)
          + 12.0*sigma.xz()*sigma.zz()*w
          + u*pow4(w)
          + 4.0*sigma.xz()*pow3(w)
        );
}

momentFuncHeader(0,1,4)
{
    moments(0,1,4)[celli] =
        m0
       *(
            3.0*v*sqr(sigma.zz())
          + 6.0*v*sigma.zz()*sqr(w)
          + 12.0*sigma.yz()*sigma.zz()*w
          + v*pow4(w)
          + 4.0*sigma.yz()*pow3(w)
        );
}

momentFuncHeader(5,1,0)
{
    moments(5,1,0)[celli] =
        m0
       *(
            15.0*v*u*sqr(sigma.xx())
          + 15.0*sigma.xy()*sqr(sigma.xx())
          + 10.0*v*sigma.xx()*pow3(u)
          + 30.0*sigma.xy()*sigma.xx()*sqr(u)
          + v*pow5(u)
          + 5.0*sigma.xy()*pow4(u)
        );
}

momentFuncHeader(5,0,1)
{
    moments(5,0,1)[celli] =
        m0
       *(
            15.0*w*u*sqr(sigma.xx())
          + 15.0*sigma.xz()*sqr(sigma.xx())
          + 10.0*w*sigma.xx()*pow3(u)
          + 30.0*sigma.xz()*sigma.xx()*sqr(u)
          + w*pow5(u)
          + 5.0*sigma.xz()*pow4(u)
        );
}

momentFuncHeader(0,5,1)
{
    moments(0,5,1)[celli] =
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

momentFuncHeader(1,5,0)
{
    moments(1,5,0)[celli] =
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

momentFuncHeader(1,0,5)
{
    moments(1,0,5)[celli] =
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

momentFuncHeader(0,1,5)
{
    moments(0,1,5)[celli] =
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

momentFuncHeader(2,2,0)
{
    moments(2,2,0)[celli] =
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

momentFuncHeader(2,0,2)
{
    moments(2,0,2)[celli] =
        m0
       *(
            2.0*sqr(sigma.xz())
          + 4.0*sigma.xz()*u*w
          + sqr(u*w)
          + sigma.zz()*sqr(u) + sigma.xx()*sqr(w)
          + sigma.xx()*sigma.zz()
        );
}

momentFuncHeader(0,2,2)
{
    moments(0,2,2)[celli] =
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
// ************************************************************************* //
