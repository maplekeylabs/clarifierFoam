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

scalar rp
(
    readScalar(settlingProperties.subDict("settlingFunction").lookup("rp"))
);

scalar rh
(
    readScalar(settlingProperties.subDict("settlingFunction").lookup("rh"))
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
