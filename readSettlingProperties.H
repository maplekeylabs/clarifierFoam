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

Info << "  Reading turbulent Schmidt number, ScT" << endl;
scalar ScT
(
    readScalar(settlingProperties.lookup("ScT"))
);
