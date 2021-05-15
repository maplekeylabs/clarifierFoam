// Calcualate the sludge density

volScalarField rho = rho_r + C*(1.0 - rho_r/rho_p);

// Create gravity vector

dimensionedVector g
(
    dimensionedVector(
        "g",
        dimAcceleration,
	9.81*g_dir
    )
);

// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevSigma(U)
  - g*(rho - rho_r)/rho
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
