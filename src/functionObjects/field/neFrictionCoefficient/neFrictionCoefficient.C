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

#include "neFrictionCoefficient.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(neFrictionCoefficient, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        neFrictionCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::neFrictionCoefficient::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Cf ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::neFrictionCoefficient::lookupOrReadVector
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

Foam::functionObjects::neFrictionCoefficient::neFrictionCoefficient
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    rhoInf_(0),
    magUInf_(0),
    wssName_("wallShearStress"),
    resultName_(name),
    patchSet_()
{
    read(dict);

    volScalarField* CfPtr
    (
        new volScalarField
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
            dimensionedScalar("0", dimless, 0.0)
        )
    );

    mesh_.objectRegistry::store(CfPtr);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::neFrictionCoefficient::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    rhoInf_     = dict.get<scalar>("rhoInf");
    magUInf_    = dict.get<scalar>("magUInf");
    wssName_    = dict.getOrDefault<word>("wallShearStress", "wallShearStress");
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


bool Foam::functionObjects::neFrictionCoefficient::execute()
{
    volScalarField& Cf =
        mesh_.lookupObjectRef<volScalarField>(resultName_);

    tmp<volVectorField> twss = lookupOrReadVector(wssName_);
    const volVectorField& wss = twss();

    if (rhoInf_ < SMALL || magUInf_ < SMALL)
    {
        WarningInFunction
            << "Freestream reference is ~0 (rhoInf = " << rhoInf_
            << ", magUInf = " << magUInf_ << "); set rhoInf and magUInf in "
            << "the function dict. " << resultName_ << " not computed." << endl;
        return false;
    }

    const scalar qInf = 0.5*rhoInf_*sqr(magUInf_);

    volScalarField::Boundary& CfBf = Cf.boundaryFieldRef();
    const volVectorField::Boundary& wssBf = wss.boundaryField();

    for (const label patchi : patchSet_)
    {
        const vectorField& Sfp   = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

        const vectorField n(-Sfp/magSfp);
        const vectorField& tau = wssBf[patchi];

        // Remove the wall-normal component, keep the tangential shear
        const vectorField tauTang(tau - n*(n & tau));

        CfBf[patchi] = mag(tauTang)/qInf;
    }

    return true;
}


bool Foam::functionObjects::neFrictionCoefficient::write()
{
    const volScalarField& Cf =
        mesh_.lookupObject<volScalarField>(resultName_);

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << Cf.name() << endl;

    Cf.write();

    const volScalarField::Boundary& CfBf = Cf.boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        const scalarField& Cfp = CfBf[patchi];

        const scalar minCf = gMin(Cfp);
        const scalar maxCf = gMax(Cfp);
        const scalar avgCf = gAverage(Cfp);

        if (Pstream::master())
        {
            Log << "    patch " << patch.name()
                << " Cf : min = " << minCf << ", max = " << maxCf
                << ", average = " << avgCf << nl;

            writeCurrentTime(file());
            file()
                << token::TAB << patch.name()
                << token::TAB << minCf
                << token::TAB << maxCf
                << token::TAB << avgCf
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //