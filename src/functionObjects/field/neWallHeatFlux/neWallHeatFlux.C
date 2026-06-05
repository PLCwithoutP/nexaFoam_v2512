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

#include "neWallHeatFlux.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(neWallHeatFlux, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        neWallHeatFlux,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::neWallHeatFlux::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Wall heat flux [W/m^2]");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::neWallHeatFlux::lookupOrReadScalar
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

Foam::functionObjects::neWallHeatFlux::neWallHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    kappaTRName_("kappaTR"),
    TTRName_("TTR"),
    Twall_(0.0),
    haveTwall_(false),
    twoTemperature_(false),
    kappaVeName_("kappaVe"),
    TVibName_("TVib"),
    resultName_(name),
    patchSet_()
{
    read(dict);

    volScalarField* qPtr
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
            dimensionedScalar("0", dimEnergy/dimTime/dimArea, 0.0)
        )
    );

    mesh_.objectRegistry::store(qPtr);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::neWallHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    kappaTRName_    = dict.getOrDefault<word>("kappaTR", "kappaTR");
    TTRName_        = dict.getOrDefault<word>("TTR", "TTR");

    // Optional solid wall temperature. If supplied, the wall-normal gradient
    // is referenced to this fixed surface temperature instead of the stored
    // boundary face value. Use this when the wall BC is a temperature-jump
    // type whose written face value collapses onto the near-wall cell under
    // -postProcess (so the gas-side gradient on disk is zero).
    haveTwall_ = dict.readIfPresent<scalar>("Twall", Twall_);

    twoTemperature_ = dict.getOrDefault<bool>("twoTemperature", false);
    kappaVeName_    = dict.getOrDefault<word>("kappaVe", "kappaVe");
    TVibName_       = dict.getOrDefault<word>("TVib", "TVib");
    resultName_     = dict.getOrDefault<word>("result", this->name());

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


bool Foam::functionObjects::neWallHeatFlux::execute()
{
    volScalarField& q =
        mesh_.lookupObjectRef<volScalarField>(resultName_);

    tmp<volScalarField> tkTR  = lookupOrReadScalar(kappaTRName_);
    tmp<volScalarField> tTTR  = lookupOrReadScalar(TTRName_);
    const volScalarField& kTR = tkTR();
    const volScalarField& TTR = tTTR();

    // Geometric 1/delta on each boundary face (cell-centre to face).
    const surfaceScalarField::Boundary& dcBf =
        mesh_.deltaCoeffs().boundaryField();

    volScalarField::Boundary& qBf = q.boundaryFieldRef();

    for (const label patchi : patchSet_)
    {
        // Near-wall cell conductivity: the thermo does not reliably populate
        // the boundary value of kappa, but the adjacent cell value is valid.
        const scalarField kw(kTR.boundaryField()[patchi].patchInternalField());

        // Explicit wall-normal gradient from the stored gas-side wall value
        // and the adjacent cell value. This avoids the patch BC's snGrad(),
        // which a jump/slip BC does not evaluate correctly under -postProcess
        // (its mixed-BC coefficients are not in the field file).
        const fvPatchScalarField& TTRp = TTR.boundaryField()[patchi];
        const scalarField Tface
        (
            haveTwall_
          ? scalarField(TTRp.size(), Twall_)  // fixed solid wall temperature
          : scalarField(TTRp)                 // stored gas wall-face value
        );
        const scalarField Tcell(TTRp.patchInternalField()); // adjacent cell T
        const scalarField& dc = dcBf[patchi];

        // Heat flux into the wall, positive = surface heating:
        //   q = -kappa * snGrad_out(T) = kappa * (T_cell - T_face) / delta
        qBf[patchi] = kw*(Tcell - Tface)*dc;
    }

    if (twoTemperature_)
    {
        tmp<volScalarField> tkVe  = lookupOrReadScalar(kappaVeName_);
        tmp<volScalarField> tTVib = lookupOrReadScalar(TVibName_);
        const volScalarField& kVe  = tkVe();
        const volScalarField& TVib = tTVib();

        for (const label patchi : patchSet_)
        {
            const scalarField kw(kVe.boundaryField()[patchi].patchInternalField());

            const fvPatchScalarField& TVp = TVib.boundaryField()[patchi];
            const scalarField Tface(TVp);
            const scalarField Tcell(TVp.patchInternalField());
            const scalarField& dc = dcBf[patchi];

            qBf[patchi] = qBf[patchi] + kw*(Tcell - Tface)*dc;
        }
    }

    return true;
}


bool Foam::functionObjects::neWallHeatFlux::write()
{
    const volScalarField& q =
        mesh_.lookupObject<volScalarField>(resultName_);

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << q.name() << endl;

    q.write();

    const volScalarField::Boundary& qBf = q.boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        const scalarField& qp = qBf[patchi];

        const scalar minQ = gMin(qp);
        const scalar maxQ = gMax(qp);
        const scalar avgQ = gAverage(qp);

        if (Pstream::master())
        {
            Log << "    patch " << patch.name()
                << " q_w : min = " << minQ << ", max = " << maxQ
                << ", average = " << avgQ << nl;

            writeCurrentTime(file());
            file()
                << token::TAB << patch.name()
                << token::TAB << minQ
                << token::TAB << maxQ
                << token::TAB << avgQ
                << endl;
        }
    }

    return true;
}


// ************************************************************************* //