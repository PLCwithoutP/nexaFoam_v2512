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


bool Foam::functionObjects::neWallHeatFlux::fieldAvailable
(
    const word& fieldName
) const
{
    // Already constructed (runtime): present in the registry.
    if (mesh_.foundObject<volScalarField>(fieldName))
    {
        return true;
    }

    // Otherwise (post-process): non-fatal check that a readable field file of
    // the right type exists in the current time directory.
    IOobject io
    (
        fieldName,
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    return io.typeHeaderOk<volScalarField>(true);
}


void Foam::functionObjects::neWallHeatFlux::addConductionTerm
(
    volScalarField::Boundary& qBf,
    const volScalarField& kappa,
    const volScalarField& T,
    const bool useWallT,
    const scalar wallT
) const
{
    const surfaceScalarField::Boundary& dcBf =
        mesh_.deltaCoeffs().boundaryField();

    for (const label patchi : patchSet_)
    {
        // Near-wall cell conductivity: the thermo does not reliably populate
        // the boundary value of kappa, but the adjacent cell value is valid.
        const scalarField kw(kappa.boundaryField()[patchi].patchInternalField());

        // Explicit wall-normal gradient from stored values, independent of the
        // patch BC's snGrad() (which a jump/slip BC does not evaluate correctly
        // under -postProcess). When a fixed wall temperature is supplied it is
        // used as the surface reference instead of the stored face value.
        const fvPatchScalarField& Tp = T.boundaryField()[patchi];
        const scalarField Tface
        (
            useWallT
          ? scalarField(Tp.size(), wallT)   // fixed solid wall temperature
          : scalarField(Tp)                  // stored wall-face value
        );
        const scalarField Tcell(Tp.patchInternalField()); // adjacent cell T
        const scalarField& dc = dcBf[patchi];

        // Accumulate; positive = heat into the wall (surface heating):
        //   q = -kappa * snGrad_out(T) = kappa * (T_cell - T_face) / delta
        qBf[patchi] += kw*(Tcell - Tface)*dc;
    }
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
    autoVib_(true),
    kappaVeName_("kappaVe"),
    TVibName_("TVib"),
    TVibWall_(0.0),
    haveTVibWall_(false),
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

    kappaVeName_ = dict.getOrDefault<word>("kappaVe", "kappaVe");
    TVibName_    = dict.getOrDefault<word>("TVib", "TVib");

    // Optional fixed wall vibrational temperature, analogous to Twall, for
    // jump/slip TVib walls under -postProcess.
    haveTVibWall_ = dict.readIfPresent<scalar>("TVibWall", TVibWall_);

    // Vibrational-conduction control:
    //   - "twoTemperature" present  -> honour it exactly (true/false)
    //   - "twoTemperature" absent   -> auto: add the vib term iff both the
    //                                  vib conductivity and vib temperature
    //                                  fields can be found at run time.
    // This lets one dictionary serve 1T, 2T, pure-gas, mixture, reacting and
    // non-reacting cases without edits.
    if (dict.found("twoTemperature"))
    {
        twoTemperature_ = dict.get<bool>("twoTemperature");
        autoVib_ = false;
    }
    else
    {
        twoTemperature_ = false;
        autoVib_ = true;
    }

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


bool Foam::functionObjects::neWallHeatFlux::execute()
{
    volScalarField& q =
        mesh_.lookupObjectRef<volScalarField>(resultName_);

    volScalarField::Boundary& qBf = q.boundaryFieldRef();

    // Reset wall patches before accumulating the conduction terms.
    for (const label patchi : patchSet_)
    {
        qBf[patchi] = 0.0;
    }

    // --- Translational-rotational conduction (always present) ---
    {
        tmp<volScalarField> tkTR = lookupOrReadScalar(kappaTRName_);
        tmp<volScalarField> tTTR = lookupOrReadScalar(TTRName_);
        addConductionTerm(qBf, tkTR(), tTTR(), haveTwall_, Twall_);
    }

    // --- Vibrational-electronic conduction (two-temperature only) ---
    // Included when explicitly requested, or (auto mode) whenever both the
    // vibrational conductivity and vibrational temperature fields exist. In
    // 1T runs the fields are absent (or published as zero), so the term is
    // skipped or contributes nothing - the object is safe in every mode.
    const bool wantVib = autoVib_ ? true : twoTemperature_;

    if (wantVib)
    {
        if (fieldAvailable(kappaVeName_) && fieldAvailable(TVibName_))
        {
            tmp<volScalarField> tkVe  = lookupOrReadScalar(kappaVeName_);
            tmp<volScalarField> tTVib = lookupOrReadScalar(TVibName_);
            addConductionTerm(qBf, tkVe(), tTVib(), haveTVibWall_, TVibWall_);
        }
        else if (!autoVib_)
        {
            WarningInFunction
                << "twoTemperature was requested but '" << kappaVeName_
                << "' and/or '" << TVibName_ << "' are unavailable; "
                << "the vibrational term is omitted." << endl;
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