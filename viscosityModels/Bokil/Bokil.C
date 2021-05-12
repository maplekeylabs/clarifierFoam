/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Bokil.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(Bokil, 0);
    addToRunTimeSelectionTable(viscosityModel, Bokil, dictionary);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Bokil::calcNu() const
{
    return pos(Cmin_-C_)*a1_*exp(b1_*C_) + pos(C_-Cmin_)*a2_*exp(b2_*C_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::Bokil::Bokil
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    BokilCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    a1_("a1", dimViscosity, BokilCoeffs_),
    b1_("b1", dimVol/dimMass, BokilCoeffs_),
    a2_("a2", dimViscosity, BokilCoeffs_),
    b2_("b2", dimVol/dimMass, BokilCoeffs_),
    Cmin_("Cmin", dimMass/dimVol, BokilCoeffs_),
    C_(U.db().lookupObject<volScalarField>("C")),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::Bokil::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    BokilCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    BokilCoeffs_.lookup("a1") >> a1_;
    BokilCoeffs_.lookup("b1") >> b1_;
    BokilCoeffs_.lookup("a2") >> a2_;
    BokilCoeffs_.lookup("b2") >> b2_;
    BokilCoeffs_.lookup("Cmin") >> Cmin_;

    return true;
}


// ************************************************************************* //
