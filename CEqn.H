volScalarField Vs = V0*(exp(-rp*(C-Cmin)) + exp(-rh*(C-Cmin)));
Vs.max(VSMALL);

fvScalarMatrix CEqn
(
    fvm::ddt(C)
    + fvm::div(phi, C)
    ==
    fvm::laplacian(turbulence->nuEff()/ScT, C)
    - Vs*(fvc::grad(C) & g)
);

CEqn.relax();
CEqn.solve();
