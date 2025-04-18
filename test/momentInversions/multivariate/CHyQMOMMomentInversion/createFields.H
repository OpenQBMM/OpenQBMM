Info<< "Reading quadratureProperties\n" << endl;

dictionary quadratureProperties(IFstream("quadratureProperties")());

labelListList nodeIndexes = quadratureProperties.lookup("nodes");
labelListList momentOrders = quadratureProperties.lookup("moments");
labelList velocityIndexes;
forAll(momentOrders[0], cmpt)
{
    velocityIndexes.append(cmpt);
}

label nNodes = nodeIndexes.size();
label nMoments = momentOrders.size();
label nDims = momentOrders[0].size();
