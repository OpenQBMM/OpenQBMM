#include "compressibleCourantNo.H"

{
    scalar realizableCo = mixing().realizableCo();

    Info<< "Max Courant Number for realizable moment transport  = "
        << realizableCo << endl;

     maxCo = min(maxCo, realizableCo);
}
