Info<< "Reading settlingProperties" << endl;

IOdictionary settlingProperties
(
    IOobject
    (
        "settlingProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< " Reading solution type" << endl;

bool includeMomentum
(
    settlingProperties.lookupOrDefault<Switch>("includeMomentum", true)
);

Info<< " Reading sludge model paremters" << endl;

dimensionedScalar rho_r
(
    dimensionedScalar(
        "rho_r",
	dimMass/dimVol,
        readScalar(settlingProperties.subDict("sludgeModel").lookup("rho_r"))
    )
);

dimensionedScalar rho_p
(
    dimensionedScalar(
        "rho_p",
        dimMass/dimVol,
        readScalar(settlingProperties.subDict("sludgeModel").lookup("rho_p"))
    )
);

Info<< "  Reading settling function parameters" << endl;

vector g_dir
(
    settlingProperties.subDict("settlingFunction").lookup("g_dir")
);

dimensionedScalar V0
(
    dimensionedScalar(
        "V0",
	dimLength/dimTime,
        readScalar(settlingProperties.subDict("settlingFunction").lookup("V0"))
    )
);

dimensionedScalar rp
(
    dimensionedScalar(
        "rp",
	dimVol/dimMass,
        readScalar(settlingProperties.subDict("settlingFunction").lookup("rp"))
    )
);

dimensionedScalar rh
(
    dimensionedScalar(
        "rh",
	dimVol/dimMass,
        readScalar(settlingProperties.subDict("settlingFunction").lookup("rh"))
    )
);

dimensionedScalar Cmin
(
    dimensionedScalar(
        "Cmin",
	dimMass/dimVol,
        readScalar(settlingProperties.subDict("settlingFunction").lookup("Cmin"))
     )
);

Info<< "  Reading scraper information" << endl;

Switch includeScrapers
(
    settlingProperties.subDict("scrapers").lookupOrDefault<Switch>("includeScrapers", false)
);

DynamicList<word, 1, 0> scraperRegions;
scalar scraperAngle = 0.0;
scalar omegaScraper = 0.0;
// Added for sloped sludge_flux zone
scalar initialScraperHeight = 0.0;
scalar finalScraperHeight = 0.0;
scalar RASDiameter = 0.0;
scalar clarifierDiameter = 0.0;

if (includeScrapers)
{
    settlingProperties.subDict("scrapers").lookup("scraperRegions") >> scraperRegions;
    omegaScraper = readScalar(settlingProperties.subDict("scrapers").lookup("omegaScraper", 0.0));
    scraperAngle = readScalar(settlingProperties.subDict("scrapers").lookup("scraperAngle"));

    initialScraperHeight = readScalar(settlingProperties.subDict("scrapers").lookup("initialScraperHeight"));
    finalScraperHeight = readScalar(settlingProperties.subDict("scrapers").lookup("finalScraperHeight"));
    RASDiameter = readScalar(settlingProperties.subDict("scrapers").lookup("RASDiameter"));
    clarifierDiameter = readScalar(settlingProperties.subDict("scrapers").lookup("clarifierDiameter"));
}

Info<< "  Reading turbulent Schmidt number, ScT" << endl;

scalar ScT
(
    readScalar(settlingProperties.lookup("ScT"))
);
