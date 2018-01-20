// interpolate weights and velocities
forAll(weights, nodei)
{
    wOwn[nodei] = fvc::interpolate(weights[nodei], own, "reconstruct(weight)");
    wNei[nodei] = fvc::interpolate(weights[nodei], nei, "reconstruct(weight)");
    UOwn[nodei] = fvc::interpolate(U[nodei], own, "reconstruct(U)");
    UNei[nodei] = fvc::interpolate(U[nodei], nei, "reconstruct(U)");
}
forAll(fluxes, mi)
{
    fluxes[mi] =
        dimensionedScalar
        (
            "0",
            moments[mi].dimensions()/dimTime,
            0.0
        );
}

dimensionedScalar zeroPhi("zeroPhi", dimVolume/dimTime, 0.0);

forAll(indicies, nodei)
{
//     surfaceScalarField phi(fvc::flux(U[nodei]));
    surfaceScalarField phiOwn(UOwn[nodei] & mesh.Sf());
    surfaceScalarField phiNei(UNei[nodei] & mesh.Sf());

    forAll(mIndicies, mi)
    {
        label oi = mIndicies[mi][0];
        label oj = mIndicies[mi][1];

        surfaceScalarField cmptOwn
        (
            wOwn[nodei]
           *pow(UOwn[nodei].component(0), oi)
           *pow(UOwn[nodei].component(1), oj)
        );

        surfaceScalarField cmptNei
        (
            wNei[nodei]
           *pow(UNei[nodei].component(0), oi)
           *pow(UNei[nodei].component(1), oj)
        );

        fluxes[mi] +=
            fvc::surfaceIntegrate
            (
                cmptOwn*max(phiOwn, zeroPhi)
              + cmptNei*min(phiNei, zeroPhi)
            );
    }
}

// forAll(moments[0].boundaryField(), patchi)
// {
//     if
//     (
//        !(
//             isA<wallFvPatch>(mesh.boundary()[patchi])
//          || isA<emptyFvPatch>(mesh.boundary()[patchi])
//         )
//     )
//     {
//         forAll(fluxes, mi)
//         {
//             fluxes[mi].boundaryFieldRef() *= -1.0;
//         }
//     }
// }
