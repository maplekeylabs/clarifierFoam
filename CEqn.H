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

if (includeScrapers)
{
/*	
    // Calculate volume of scraper regions
    forAll(scraperRegions, iRegion)
    {
        const label zoneID = mesh.cellZones().findZoneID(scraperRegions[iRegion]);

        if (zoneID == -1)
        {
            break;
        }

        const labelList& cells = mesh.cellZones()[zoneID];

        forAll(cells, celli)
        {
            scraperVolume += mesh.V()[cells[celli]];
        }
    }
    //reduce(scraperVolume, sumOp<scalar>());
    Info<< "  Scaper region volume: " << scraperVolume << endl;

    // Calculate scraper flow rate
    scraperFlow = scraperVolume*omegaScraper/2.0/constant::mathematical::pi;
    Info<< "  Scraper volume flow rate: " << scraperFlow << endl;
*/
    // Calculate the scraper velocity field
    forAll(scraperRegions, iRegion)
    {
        const label zoneID = mesh.cellZones().findZoneID(scraperRegions[iRegion]);

        if (zoneID == -1)
        {
            break;
        }

              const labelList& cells = mesh.cellZones()[zoneID];

        forAll(cells, celli)
        {
	    // Divide diameters by 2 to avoid confusing code:
	    scalar r_RAS = RASDiameter / 2;
	    scalar r_clarifier = clarifierDiameter / 2;
	    
	    // Find coordinated control volume cell center
            label i = cells[celli];
            vector pos = mesh.C()[i];

	    // Added RAS outlet radius offset
            scalar r = Foam::sqrt(Foam::pow(pos[0], 2.0) + Foam::pow(pos[1], 2.0));
            scalar scraperMag = r/Foam::cos(scraperAngle);
            vector scraperDir(-pos[0]/scraperMag, -pos[1]/scraperMag, -r*Foam::tan(scraperAngle)/scraperMag);

	    // Calculate the local height of the clarifier floor
	    scalar r_bar = r - r_RAS;
            scalar h_floor = r_bar*Foam::tan(scraperAngle);

            // Calculate the height of the sludge_flux zone
	    scalar h_scraper = initialScraperHeight + (r_bar / (r_clarifier - r_RAS) ) * (finalScraperHeight - initialScraperHeight);
	    
	    // Calculate the local scraper height 
	    scalar h = h_scraper - h_floor;

	    //Calculate the scraper velocity 
            Vscraper[i] = scraperFlow/2.0/constant::mathematical::pi/(r_bar + r_RAS)/h*scraperDir;
        }
    }
}

Info << "Scraper volume: " << scraperVolume << endl;
Info << "Scraper volume flow rate: " << scraperFlow << endl;

// Calculate phi field for scraper velocity (phisc)
surfaceScalarField phisc
(
    "phisc",
    fvc::interpolate(Vscraper) & mesh.Sf()
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
    + fvOptions(C)
);

CEqn.relax();
fvOptions.contrain(CEqn);
CEqn.solve();
