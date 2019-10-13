    Info<< "Reading quadratureProperties\n" << endl;

    dictionary quadratureProperties(IFstream("quadratureProperties")());

    labelListList nodeIndexes = quadratureProperties.lookup("nodes");
    labelListList momentOrders = quadratureProperties.lookup("moments");

    labelList velocityIndexes;
    
    for (label i = 1; i < momentOrders[0].size(); i++)
    {
        velocityIndexes.append(i);
    }

    label nNodes = nodeIndexes.size();
    label nMoments = momentOrders.size();
    label nDims = momentOrders[0].size();
