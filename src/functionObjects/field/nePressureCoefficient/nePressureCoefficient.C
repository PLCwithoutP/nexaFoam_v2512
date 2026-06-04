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

#include "nePressureCoefficient.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(nePressureCoefficient, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        nePressureCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::nePressureCoefficient::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Cp ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::nePressureCoefficient::lookupOrReadScalar
(
    const word& fieldName
) const
{
    // Runtime path: the field is live in the registry. Wrap it in a
    // non-owning tmp so we do not copy or delete it.
    if (mesh_.foundObject<volScalarField>(fieldName))
    {
        return tmp<volScalarField>
        (
            mesh_.lookupObject<volScalarField>(fieldName)
        );
    }

    // Offline (post-process) path: build an owning tmp by reading the field
    // from the current time directory. NO_REGISTER keeps it out of the
    // registry so it is freed when this tmp goes out of scope.
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::nePressureCoefficient::nePressureCoefficient
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
    resultName_(name),
    patchSet_()
{
    read(dict);

    // Create the result field and hand ownership to the registry so it
    // persists across time steps and can be found by execute(), write()
    // and any downstream sampler.
    volScalarField* CpPtr
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

    mesh_.objectRegistry::store(CpPtr);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::nePressureCoefficient::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    rhoInf_     = dict.get<scalar>("rhoInf");
    magUInf_    = dict.get<scalar>("magUInf");
    resultName_ = dict.getOrDefault<word>("result", this->name());

    // Resolve the wall patches to process.
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    patchSet_ = pbm.patchSet(dict.getOrDefault<wordRes>("patches", wordRes()));

    if (patchSet_.empty())
    {
        // No 'patches' entry: process every wall patch.
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
        // 'patches' given: keep only wall patches, warn about the rest.
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


bool Foam::functionObjects::nePressureCoefficient::execute()
{
    // Mutable handle to the registered result field (no const_cast needed).
    volScalarField& Cp =
        mesh_.lookupObjectRef<volScalarField>(resultName_);

    // Wall pressure only; the freestream reference is user-specified.
    tmp<volScalarField> tp = lookupOrReadScalar("p");
    const volScalarField& p = tp();

    if (rhoInf_ < SMALL || magUInf_ < SMALL)
    {
        WarningInFunction
            << "Freestream reference is ~0 (rhoInf = " << rhoInf_
            << ", magUInf = " << magUInf_ << "); set rhoInf and magUInf in "
            << "the function dict. " << resultName_ << " not computed." << endl;
        return false;
    }

    const scalar qInf = 0.5*rhoInf_*sqr(magUInf_);

    volScalarField::Boundary& CpBf = Cp.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p.boundaryField();

    for (const label patchi : patchSet_)
    {
        CpBf[patchi] = pBf[patchi]/qInf;
    }

    return true;
}


bool Foam::functionObjects::nePressureCoefficient::write()
{
    const volScalarField& Cp =
        mesh_.lookupObject<volScalarField>(resultName_);

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << Cp.name() << endl;

    Cp.write();

    const volScalarField::Boundary& CpBf = Cp.boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        const scalarField& Cpp = CpBf[patchi];

        const scalar minCp = gMin(Cpp);
        const scalar maxCp = gMax(Cpp);
        const scalar avgCp = gAverage(Cpp);

        if (Pstream::master())
        {
            Log << "    patch " << patch.name()
                << " Cp : min = " << minCp << ", max = " << maxCp
                << ", average = " << avgCp << nl;

            writeCurrentTime(file());
            file()
                << token::TAB << patch.name()
                << token::TAB << minCp
                << token::TAB << maxCp
                << token::TAB << avgCp
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //