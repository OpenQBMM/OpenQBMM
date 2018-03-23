#include "compressibleCourantNo.H"

{
    scalar realizableCo = populationBalance().realizableCo();

    Info<< "Max Courant Number for realizable moment transport  = "
        << realizableCo << endl;

    maxCo = min(maxCo, realizableCo);
}
