Thetap = dimensionedScalar("zero", sqr(dimVelocity), 0.0);

forAll(momentOrders[0], cmpti)
{
    labelList order1(nDims, 0);
    labelList order2(nDims, 0);
    order1[cmpti] = 1;
    order2[cmpti] = 2;

    volScalarField meanU("meanU", moments(order1)/max(alpha, residualAlpha));
    Up.replace(cmpti, meanU);

    Thetap +=
        max
        (
            moments(order2)/max(alpha, residualAlpha) - sqr(meanU),
            dimensionedScalar("zero", sqr(dimVelocity), 0.0)
        );
}
