fvScalarMatrix CEqn
(
    fvm::ddt(C)
    + fvm::div(phi, C)
    ==
    fvm::laplacian(turbulence->nuEff()/ScT, C)
);

CEqn.relax();
CEqn.solve();
