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

Info<< " Reading settling function parameters" << endl;

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

Info<< "  Reading turbulent Schmidt number, ScT" << endl;

scalar ScT
(
    readScalar(settlingProperties.lookup("ScT"))
);
