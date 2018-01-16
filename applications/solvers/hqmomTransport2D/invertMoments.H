const volScalarField& m0 = moments[0];
forAll(m0, celli)
{
    multivariateMomentSet m(10, mIndicies, "R");

    forAll(m, mi)
    {
        m[mi] = moments[mi][celli];
    }
//     Info<<"cell " << celli << ": " << m << ", " << endl;
    momentInverter.invert(m);
    const mappedList<scalar>& w = momentInverter.weights();
    const mappedList<vector>& u = momentInverter.abscissae();

//     Info<< w << nl << u << endl;

    forAll(w, nodei)
    {
        weights[nodei][celli] = w[nodei];
        U[nodei][celli] = u[nodei];
    }
}

const volScalarField::Boundary& bf = moments[0].boundaryField();
forAll(bf, patchi)
{
    const fvPatchScalarField& m0Patch = bf[patchi];
    forAll(m0Patch, facei)
    {
        multivariateMomentSet m(10, mIndicies, "R");

        forAll(m, mi)
        {
            m[mi] = moments[mi].boundaryField()[patchi][facei];
        }

        momentInverter.invert(m);
        const mappedList<scalar>& w = momentInverter.weights();
        const mappedList<vector>& u = momentInverter.abscissae();

        forAll(w, nodei)
        {
            weights[nodei].boundaryFieldRef()[patchi][facei] = w[nodei];
            U[nodei].boundaryFieldRef()[patchi][facei] = u[nodei];
        }
    }
}
