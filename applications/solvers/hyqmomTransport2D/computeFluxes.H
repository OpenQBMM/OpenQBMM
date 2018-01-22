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
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& currPatch = mesh.boundary()[patchi];
        if (isA<wallFvPatch>(currPatch))
        {
            const vectorField& bfSf = mesh.Sf().boundaryField()[patchi];
            vectorField bfNorm = bfSf/mag(bfSf);

            scalarField& bfwOwn = wOwn[nodei].boundaryFieldRef()[patchi];
            scalarField& bfwNei = wNei[nodei].boundaryFieldRef()[patchi];
            vectorField& bfUOwn = UOwn[nodei].boundaryFieldRef()[patchi];
            vectorField& bfUNei = UNei[nodei].boundaryFieldRef()[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                bfwOwn[facei] = weights[nodei][faceCelli];
                bfUOwn[facei] = U[nodei][faceCelli];

                bfwNei[facei] = bfwOwn[facei];
                bfUNei[facei] = bfUOwn[facei]
                  - 2.0*(bfUOwn[facei] & bfNorm[facei])
                   *bfNorm[facei];
            }
        }
    }

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

// // AP: Fix for wall interpolation (assuming reflective walls)
//   // This has to be generalized. At the moment it is first-order accurate
//   // and it hard-codes the reflective BC!
//   const fvPatchList& patches = mesh.boundary();
//
//   forAll(patches, patchi)
//   {
//       const fvPatch& currPatch = patches[patchi];
//
//       if (isA<wallFvPatch😠currPatch))
//       {
//     scalarField& wOwnw = weightOwn[nI].boundaryField()[patchi];
//     scalarField& wNeiw = weightNei[nI].boundaryField()[patchi];
//
//     vectorField& aOwnw = abscissaOwn[nI].boundaryField()[patchi];
//     vectorField& aNeiw = abscissaNei[nI].boundaryField()[patchi];
//
//     vectorField faceNormalVectors =
//       mesh.Sf().boundaryField()[patchi]
//         /mesh.magSf().boundaryField()[patchi];
//
//     forAll(currPatch, faceI)
//     {
//         label faceCelli = currPatch.faceCells()[faceI];
//
//         wOwnw[faceI] = weight[faceCelli];
//         aOwnw[faceI] = abscissa[faceCelli];
//
//         wNeiw[faceI] = wOwnw[faceI];
//
//         aNeiw[faceI] = aOwnw[faceI]
//       - 2.0*(aOwnw[faceI] & faceNormalVectors[faceI])
//       *faceNormalVectors[faceI];
//     }
//       }
//   }
//     }
