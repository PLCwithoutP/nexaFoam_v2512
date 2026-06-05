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

#include "neStantonNumber.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(neStantonNumber, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        neStantonNumber,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::neStantonNumber::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "St ()");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::neStantonNumber::lookupOrReadScalar
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::neStantonNumber::neStantonNumber
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
    qName_("wallHeatFlux"),
    resultName_(name),
    patchSet_()
{
    read(dict);

    volScalarField* StPtr
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

    mesh_.objectRegistry::store(StPtr);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::neStantonNumber::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    rhoInf_     = dict.get<scalar>("rhoInf");
    magUInf_    = dict.get<scalar>("magUInf");
    qName_      = dict.getOrDefault<word>("wallHeatFlux", "wallHeatFlux");
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


bool Foam::functionObjects::neStantonNumber::execute()
{
    volScalarField& St =
        mesh_.lookupObjectRef<volScalarField>(resultName_);

    tmp<volScalarField> tq = lookupOrReadScalar(qName_);
    const volScalarField& q = tq();

    if (rhoInf_ < SMALL || magUInf_ < SMALL)
    {
        WarningInFunction
            << "Freestream reference is ~0 (rhoInf = " << rhoInf_
            << ", magUInf = " << magUInf_ << "); set rhoInf and magUInf in "
            << "the function dict. " << resultName_ << " not computed." << endl;
        return false;
    }

    const scalar qRef = 0.5*rhoInf_*pow3(magUInf_);

    volScalarField::Boundary& StBf = St.boundaryFieldRef();
    const volScalarField::Boundary& qBf = q.boundaryField();

    for (const label patchi : patchSet_)
    {
        StBf[patchi] = qBf[patchi]/qRef;
    }

    return true;
}


bool Foam::functionObjects::neStantonNumber::write()
{
    const volScalarField& St =
        mesh_.lookupObject<volScalarField>(resultName_);

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << St.name() << endl;

    St.write();

    const volScalarField::Boundary& StBf = St.boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        const scalarField& Stp = StBf[patchi];

        const scalar minSt = gMin(Stp);
        const scalar maxSt = gMax(Stp);
        const scalar avgSt = gAverage(Stp);

        if (Pstream::master())
        {
            Log << "    patch " << patch.name()
                << " St : min = " << minSt << ", max = " << maxSt
                << ", average = " << avgSt << nl;

            writeCurrentTime(file());
            file()
                << token::TAB << patch.name()
                << token::TAB << minSt
                << token::TAB << maxSt
                << token::TAB << avgSt
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //