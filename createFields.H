#include "createRDeltaT.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field C\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Vs\n" << endl;
volVectorField Vs
(
    IOobject
    (
        "Vs",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rho", dimDensity, scalar(0.0))
);

Info<< "Creating field Vscraper\n" << endl;
volVectorField Vscraper
(
    IOobject
    (
        "Vscraper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Vscraper", dimVelocity, vector(0.0, 0.0, 0.0))
);

// Calculate the scraper flux parameters
scalar scraperVolume = 0.0;
scalar scraperFlow = 0.0;

if (includeScrapers)
{
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
    reduce(scraperVolume, sumOp<scalar>());
    Info<< "  Scaper region volume: " << scraperVolume << endl;

    // Calculate scraper flow rate
    scraperFlow = scraperVolume*omegaScraper/2.0/constant::mathematical::pi;
    Info<< "  Scraper volume flow rate: " << scraperFlow << endl;

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
            label i = cells[celli];
            vector pos = mesh.C()[i];
            scalar r = Foam::sqrt(Foam::pow(pos[0], 2.0) + Foam::pow(pos[1], 2.0));
	          scalar scraperMag = r/Foam::cos(scraperAngle);
            vector scraperDir(-pos[0]/scraperMag, -pos[1]/scraperMag, -r*Foam::tan(scraperAngle)/scraperMag);
            Vscraper[i] = scraperFlow/2.0/constant::mathematical::pi/r/scraperHeight*scraperDir;
        }
    }
}

// Calculate phi field for scraper velocity (phisc)
surfaceScalarField phisc
(
    "phisc",
    fvc::interpolate(Vscraper) & mesh.Sf()
);

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
