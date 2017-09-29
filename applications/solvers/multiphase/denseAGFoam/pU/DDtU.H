{
    DDtU1 =
        fvc::ddt(U1)
      + fvc::div(h2f*phi1, U1, "div(" + phi1.name() + "," + U1.name() + ")")
      - fvc::div(h2f*phi1)*U1;

    DDtU2 =
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2;
}
