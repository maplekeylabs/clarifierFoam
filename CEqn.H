// Limit the concentration field to be positive
C.max(0.0);

// Calculate the settling velocity
Vs = V0*(exp(-rp*(C-Cmin)) + exp(-rh*(C-Cmin)))*g_dir;

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
    ==
    fvm::laplacian(turbulence->nuEff() + (1.0/ScT - 1.0)*turbulence->nut(), C)
);

CEqn.relax();
CEqn.solve();
