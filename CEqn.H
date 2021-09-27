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

// Calculate the scraper flux velocity
if (includeScrapers)
{
    forAll(scraperRegions, iRegion)
    {
        const label zoneID = mesh.cellZones().findZoneID(scraperRegions[iRegion]);

        if (zoneID == -1)
        {
            FatalErrorIn("CEqn.H")
                << "Cannot find cellZone " << scraperRegions[iRegion] << endl
		<< exit(FatalError);
        }

	const labelList& cells = mesh.cellZones()[zoneID];

        forAll(cells, celli)
        {
            label i = cells[celli];
            vector pos = mesh.C()[i];
            scalar r = Foam::sqrt(Foam::pow(pos[0], 2.0) + Foam::pow(pos[1], 2.0));
            vector scraperDir(-pos[0], pos[1], r*Foam::tan(scraperAngle));
	    scraperDir /= Foam::mag(scraperDir);
            Vscraper[i] = scraperFlow/2.0/constant::mathematical::pi/r/scraperHeight*scraperDir;
        }
    }
}

// Define equation for concentration
fvScalarMatrix CEqn
(
    fvm::ddt(C)
    + fvm::div(phi, C)
    + fvm::div(phis, C)
    ==
    fvm::laplacian(turbulence->nut()/ScT, C)
);

CEqn.relax();
CEqn.solve();
