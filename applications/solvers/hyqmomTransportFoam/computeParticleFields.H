{
    const volScalarField& m0 = moments(0);
    Thetap = dimensionedScalar("zero", sqr(dimVelocity), 0.0);

    forAll(velocityIndexes, cmpti)
    {
        labelList order1(nDims, 0);
        labelList order2(nDims, 0);
        order1[velocityIndexes[cmpti]] = 1;
        order2[velocityIndexes[cmpti]] = 2;

        volScalarField meanU("meanU", moments(order1)/max(m0, 1e-10));
        Up.replace(cmpti, meanU);

        Thetap +=
            max
            (
                moments(order2)/max(m0, 1e-10) - sqr(meanU),
                dimensionedScalar("zero", sqr(dimVelocity), 0.0)
            );
    }
    Thetap /= 3.0;
}
