/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 Alberto Passalacqua
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class weightType, class abscissaType, class sigmaType>
const Foam::word&
Foam::quadratureNode<weightType, abscissaType, sigmaType>::name() const
{
    return name_;
}

template <class weightType, class abscissaType, class sigmaType>
inline bool Foam::quadratureNode<weightType, abscissaType, sigmaType>::
extended() const
{
    return extended_;
}

template <class weightType, class abscissaType, class sigmaType>
Foam::label
Foam::quadratureNode<weightType, abscissaType, sigmaType>::
nSecondaryNodes() const
{
    return nSecondaryNodes_;
}

template <class weightType, class abscissaType, class sigmaType>
const weightType&
Foam::quadratureNode<weightType, abscissaType, sigmaType>::
primaryWeight() const
{
    return weight_;
}

template <class weightType, class abscissaType, class sigmaType>
weightType&
Foam::quadratureNode<weightType, abscissaType, sigmaType>::primaryWeight()
{
    return weight_;
}

template <class weightType, class abscissaType, class sigmaType>
const abscissaType&
Foam::quadratureNode<weightType, abscissaType, sigmaType>::
primaryAbscissa() const
{
    return abscissa_;
}

template <class weightType, class abscissaType, class sigmaType>
abscissaType&
Foam::quadratureNode<weightType, abscissaType, sigmaType>::primaryAbscissa()
{
    return abscissa_;
}

template <class weightType, class abscissaType, class sigmaType>
const Foam::PtrList<weightType>&
Foam::quadratureNode<weightType, abscissaType, sigmaType>::
secondaryWeights() const
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access secondary weights of a quadrature node" << nl
            << "    not of extended type. "
            << abort(FatalError);
    }

    return secondaryWeights_;
}

template <class weightType, class abscissaType, class sigmaType>
Foam::PtrList<weightType>&
Foam::quadratureNode<weightType, abscissaType, sigmaType>::
secondaryWeights()
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access secondary weights of a quadrature node" << nl
            << "    not of extended type. "
            << abort(FatalError);
    }
    return secondaryWeights_;
}

template <class weightType, class abscissaType, class sigmaType>
const Foam::PtrList<abscissaType>&
Foam::quadratureNode<weightType, abscissaType, sigmaType>
::secondaryAbscissae() const
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access secondary abscissae of a quadrature" << nl
            << "   node not of extended type. "
            << abort(FatalError);
    }

    return secondaryAbscissae_;
}

template <class weightType, class abscissaType, class sigmaType>
Foam::PtrList<abscissaType>&
Foam::quadratureNode<weightType, abscissaType, sigmaType>
::secondaryAbscissae()
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access secondary abscissae of a quadrature" << nl
            << "   node not of extended type. "
            << abort(FatalError);
    }

    return secondaryAbscissae_;
}

template <class weightType, class abscissaType, class sigmaType>
const sigmaType&
Foam::quadratureNode<weightType, abscissaType, sigmaType>
::sigma() const
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access sigma field of a quadrature node not" << nl
            << "   of extended type. "
            << abort(FatalError);
    }

    return sigma_();
}

template <class weightType, class abscissaType, class sigmaType>
sigmaType&
Foam::quadratureNode<weightType, abscissaType, sigmaType>::sigma()
{
    if (!extended_)
    {
        FatalErrorInFunction
            << "Attempt to access sigma field of a quadrature node not" << nl
            << "   of extended type. "
            << abort(FatalError);
    }

    return sigma_();
}


// ************************************************************************* //
