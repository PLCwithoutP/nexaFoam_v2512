/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

\*---------------------------------------------------------------------------*/

#include "neWallShearStress.H"
#include "volFields.H"
#include "fvc.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(neWallShearStress, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        neWallShearStress,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::neWallShearStress::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Wall shear stress magnitude [Pa]");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::neWallShearStress::lookupOrReadScalar
(
    const word& fieldName
) const
{
    if (mesh_.foundObject<volScalarField>(fieldName))
    {
        return tmp<volScalarField>
        (
            mesh_.lookupObject<volScalarField>(fieldName)
        );
    }

    return tmp<volScalarField>::New
    (
        IOobject
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_
    );
}


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::neWallShearStress::lookupOrReadVector
(
    const word& fieldName
) const
{
    if (mesh_.foundObject<volVectorField>(fieldName))
    {
        return tmp<volVectorField>
        (
            mesh_.lookupObject<volVectorField>(fieldName)
        );
    }

    return tmp<volVectorField>::New
    (
        IOobject
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::neWallShearStress::neWallShearStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    muName_("muLam"),
    UName_("U"),
    resultName_(name),
    patchSet_()
{
    read(dict);

    // Result field: wall traction vector [Pa], zero in the interior.
    volVectorField* wssPtr
    (
        new volVectorField
        (
            IOobject
            (
                resultName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector(dimPressure, Zero)
        )
    );

    mesh_.objectRegistry::store(wssPtr);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::neWallShearStress::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    muName_     = dict.getOrDefault<word>("mu", "muLam");
    UName_      = dict.getOrDefault<word>("U", "U");
    resultName_ = dict.getOrDefault<word>("result", this->name());

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    patchSet_ = pbm.patchSet(dict.getOrDefault<wordRes>("patches", wordRes()));

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }
    }
    else
    {
        labelHashSet wallSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                wallSet.insert(patchi);
            }
            else
            {
                WarningInFunction
                    << "Requested patch '" << pbm[patchi].name()
                    << "' is not a wall; skipping it." << endl;
            }
        }
        patchSet_ = wallSet;
    }

    return true;
}


bool Foam::functionObjects::neWallShearStress::execute()
{
    volVectorField& wss =
        mesh_.lookupObjectRef<volVectorField>(resultName_);

    tmp<volScalarField> tmu = lookupOrReadScalar(muName_);
    tmp<volVectorField> tU  = lookupOrReadVector(UName_);
    const volScalarField& mu = tmu();
    const volVectorField& U  = tU();

    // Laminar deviatoric viscous stress (= laminar devRhoReff)
    const tmp<volSymmTensorField> tReff
    (
        -mu*dev(twoSymm(fvc::grad(U)))
    );
    const volSymmTensorField::Boundary& ReffBf = tReff().boundaryField();

    volVectorField::Boundary& wssBf = wss.boundaryFieldRef();

    for (const label patchi : patchSet_)
    {
        const vectorField& Sfp   = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

        // Traction on the wall using the inward unit normal -Sf/|Sf|
        wssBf[patchi] = (-Sfp/magSfp) & ReffBf[patchi];
    }

    return true;
}


bool Foam::functionObjects::neWallShearStress::write()
{
    const volVectorField& wss =
        mesh_.lookupObject<volVectorField>(resultName_);

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << wss.name() << endl;

    wss.write();

    const volVectorField::Boundary& wssBf = wss.boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        const scalarField magTau(mag(wssBf[patchi]));

        const scalar minTau = gMin(magTau);
        const scalar maxTau = gMax(magTau);
        const scalar avgTau = gAverage(magTau);

        if (Pstream::master())
        {
            Log << "    patch " << patch.name()
                << " |tau_w| : min = " << minTau << ", max = " << maxTau
                << ", average = " << avgTau << nl;

            writeCurrentTime(file());
            file()
                << token::TAB << patch.name()
                << token::TAB << minTau
                << token::TAB << maxTau
                << token::TAB << avgTau
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //