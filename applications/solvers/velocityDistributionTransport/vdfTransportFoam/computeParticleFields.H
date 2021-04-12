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
        labelList order200(orderZero);
        labelList order110(orderZero);
        labelList order101(orderZero);
        labelList order020(orderZero);
        labelList order011(orderZero);
        labelList order002(orderZero);

        order200[velocityIndexes[0]] = 2;

        Sigmap.replace
        (
            symmTensor::XX,
            moments(order200)/m0 - sqr(Up.component(0))
        );

        if (velocityIndexes.size() > 1)
        {
            order110[velocityIndexes[0]] = 1;
            order110[velocityIndexes[1]] = 1;
            order020[velocityIndexes[1]] = 2;
            Sigmap.replace
            (
                symmTensor::XY,
                moments(order110)/m0 - Up.component(0)*Up.component(1)
            );
            Sigmap.replace
            (
                symmTensor::YY,
                moments(order020)/m0 - sqr(Up.component(1))
            );
        }

        if (velocityIndexes.size() > 2)
        {
            order101[velocityIndexes[0]] = 1;
            order101[velocityIndexes[2]] = 1;
            order011[velocityIndexes[1]] = 1;
            order011[velocityIndexes[2]] = 1;
            order002[velocityIndexes[2]] = 2;

            Sigmap.replace
            (
                symmTensor::XZ,
                moments(order101)/m0 - Up.component(0)*Up.component(2)
            );

            Sigmap.replace
            (
                symmTensor::YZ,
                moments(order011)/m0 - Up.component(1)*Up.component(2)
            );

            Sigmap.replace
            (
                symmTensor::ZZ,
                moments(order002)/m0 - sqr(Up.component(2))
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
