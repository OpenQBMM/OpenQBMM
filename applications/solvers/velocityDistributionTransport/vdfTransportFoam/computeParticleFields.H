{
    const labelList& velocityIndexes = nodes[0].velocityIndexes();
    labelList orderZero(momentOrders[0].size(), 0);
    volScalarField m0(moments(0));
    m0.max(1e-10);

    forAll(velocityIndexes, cmpt)
    {
        labelList orderOne(orderZero);
        orderOne[velocityIndexes[cmpt]] = 1;

        volScalarField meanU(moments(orderOne)/m0);
        Up.replace(cmpt, meanU);
    }

    if (computeVariance)
    {
        Sigmap.replace
        (
            symmTensor::XX,
            moments(2)/m0 - sqr(Up.component(0))
        );
        if (velocityIndexes.size() > 1)
        {
            Sigmap.replace
            (
                symmTensor::XY,
                moments(1, 1)/m0 - Up.component(0)*Up.component(1)
            );
            Sigmap.replace
            (
                symmTensor::YY,
                moments(0, 2)/m0 - sqr(Up.component(1))
            );
        }
        if (velocityIndexes.size() > 2)
        {
            Sigmap.replace
            (
                symmTensor::XZ,
                moments(1, 0, 1)/m0 - Up.component(0)*Up.component(2)
            );
            Sigmap.replace
            (
                symmTensor::YZ,
                moments(0, 1, 1)/m0 - Up.component(1)*Up.component(2)
            );
            Sigmap.replace
            (
                symmTensor::ZZ,
                moments(0, 0, 2)/m0 - sqr(Up.component(2))
            );
        }
        Thetap = tr(Sigmap)/velocityIndexes.size();
    }

    if (sizeIndex != -1 && nodes[0].lengthBased())
    {
        labelList orderOne(orderZero);
        orderOne[sizeIndex] = 1;

        dMean.ref() = moments(orderOne)/m0;
    }
}
