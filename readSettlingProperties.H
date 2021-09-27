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
scalar scraperHeight = 0.0;
scalar scraperAngle = 0.0;
scalar scraperVolume = 0.0;
scalar scraperFlow = 0.0;

if (includeScrapers)
{
    settlingProperties.subDict("scrapers").lookup("scraperRegions") >> scraperRegions;

    scraperHeight = readScalar(settlingProperties.subDict("scrapers").lookup("scraperHeight"));
    scraperAngle = readScalar(settlingProperties.subDict("scrapers").lookup("scraperAngle"));

    scraperVolume = 0.0;
    forAll(scraperRegions, iRegion)
    {
        const label zoneID = mesh.cellZones().findZoneID(scraperRegions[iRegion]);

        if (zoneID == -1)
        {
            FatalErrorIn("readSettlingProperties.H")
                << "Cannot find cellZone " << scraperRegions[iRegion] << endl
		<< exit(FatalError);
        }

	const labelList& cells = mesh.cellZones()[zoneID];

        forAll(cells, celli)
        {
            scraperVolume += mesh.V()[celli];
        }
    }
    Info<< "    Scaper region volume: " << scraperVolume << endl;

    scalar omegaScraper
    (
        settlingProperties.subDict("scrapers").lookupOrDefault<scalar>("omegaScraper", 0.0)
    );

    scraperFlow = scraperVolume*omegaScraper/2.0/constant::mathematical::pi;
    Info<< "    Scraper volume flow rate: " << scraperFlow << endl;
}

Info<< "  Reading turbulent Schmidt number, ScT" << endl;

scalar ScT
(
    readScalar(settlingProperties.lookup("ScT"))
);
