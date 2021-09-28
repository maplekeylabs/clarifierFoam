// Limit the concentration field to be positive
C.max(0.0);
C.correctBoundaryConditions();

// Calculate the settling velocity
Vs = V0*(exp(-rh*(C-Cmin)) - exp(-rp*(C-Cmin)))*g_dir;

// Constrain settling velocity to be in gravity direction
forAll(Vs, celli)
{
    if ((Vs[celli] & g_dir) < 0.0)
    {
        Vs[celli] = vector(0,0,0);
    }
}

// Correct settling velocity boundary conditions
Vs.correctBoundaryConditions();

// Calculate phi field for settling velocity (phis)
surfaceScalarField phis
(
    "phis",
    fvc::interpolate(Vs) & mesh.Sf()
);

// Define equation for concentration
fvScalarMatrix CEqn
(
    fvm::ddt(C)
    + fvm::div(phi, C)
    + fvm::div(phis, C)
    + fvm::div(phisc, C)
    ==
    fvm::laplacian(turbulence->nut()/ScT, C)
);

CEqn.relax();
CEqn.solve();
