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

#include "neVibWallHeatFlux.H"
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
    defineTypeNameAndDebug(neVibWallHeatFlux, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        neVibWallHeatFlux,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::neVibWallHeatFlux::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Vibrational-electronic wall heat flux [W/m^2]");

    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::neVibWallHeatFlux::lookupOrReadScalar
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

    // kappaVe / TVib are mandatory for this object (it exists specifically
    // to report the vibrational term), so MUST_READ is intentionally left
    // to fail fatally here rather than being probed first, unlike the
    // optional-term handling in neWallHeatFlux.
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


void Foam::functionObjects::neVibWallHeatFlux::addConductionTerm
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
        // Wall-FACE conductivity. he2TThermo::kappaVib(patchi) populates every
        // boundary face per species from patch-local p/TTR/TVib/Y, so the
        // patch value is valid and is what the solver's Stage 3 laplacian
        // uses. patchInternalField() would take kappa from the adjacent CELL:
        // CvVib is exponential in T, so with TVib_cell ~ 10^3 K against a
        // 300 K wall that is a several-hundred-fold overestimate, and it makes
        // the reported flux disagree with the flux the solver applied.
        const scalarField kw(kappa.boundaryField()[patchi]);

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

Foam::functionObjects::neVibWallHeatFlux::neVibWallHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
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

bool Foam::functionObjects::neVibWallHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    kappaVeName_ = dict.getOrDefault<word>("kappaVe", "kappaVe");
    TVibName_    = dict.getOrDefault<word>("TVib", "TVib");

    // Optional solid wall vibrational temperature. If supplied, the
    // wall-normal gradient is referenced to this fixed surface temperature
    // instead of the stored boundary face value. Use this when the wall BC
    // is a temperature-jump type whose written face value collapses onto
    // the near-wall cell under -postProcess (so the gas-side gradient on
    // disk is zero).
    haveTVibWall_ = dict.readIfPresent<scalar>("TVibWall", TVibWall_);

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


bool Foam::functionObjects::neVibWallHeatFlux::execute()
{
    volScalarField& q =
        mesh_.lookupObjectRef<volScalarField>(resultName_);

    volScalarField::Boundary& qBf = q.boundaryFieldRef();

    // Reset wall patches before accumulating the conduction term.
    for (const label patchi : patchSet_)
    {
        qBf[patchi] = 0.0;
    }

    // --- Vibrational-electronic conduction (sole term) ---
    tmp<volScalarField> tkVe  = lookupOrReadScalar(kappaVeName_);
    tmp<volScalarField> tTVib = lookupOrReadScalar(TVibName_);
    addConductionTerm(qBf, tkVe(), tTVib(), haveTVibWall_, TVibWall_);

    return true;
}


bool Foam::functionObjects::neVibWallHeatFlux::write()
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
                << " q_w,vib : min = " << minQ << ", max = " << maxQ
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
