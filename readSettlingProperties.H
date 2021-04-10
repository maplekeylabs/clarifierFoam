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

vector g
(
    settlingProperties.subDict("settlingFunction").lookup("g")
);

scalar V0
(
    readScalar(settlingProperties.subDict("settlingFunction").lookup("V0"))
);

scalar rp
(
    readScalar(settlingProperties.subDict("settlingFunction").lookup("rp"))
);

scalar rh
(
    readScalar(settlingProperties.subDict("settlingFunction").lookup("rh"))
);

scalar Cmin
(
    readScalar(settlingProperties.subDict("settlingFunction").lookup("Cmin"))
);

Info<< "  Reading turbulent Schmidt number, ScT" << endl;

scalar ScT
(
    readScalar(settlingProperties.lookup("ScT"))
);
