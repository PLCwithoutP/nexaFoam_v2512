#include "createRDeltaT.H"
#include "thermochemistryCheck.H"
#include "createMixture.H"

volScalarField& p = thermo2T.p();
volScalarField& TVib = thermo2T.TVib();
volScalarField& TTR = thermo2T.TTR();
const volScalarField& psi = thermo2T.psi();

if (!inviscid)
{
    Info<< "Viscous analysis detected.\n" << endl;
}
else
{
    Info<< "Inviscid analysis detected.\n" << endl;
}

PtrList<volScalarField>& Y_species = mix.Y();

const bool solveSpeciesVib = mixtureCheck && use2T;
const bool includeAtomicVe = thermo2T.hasElectronicVibCoupling();

List<bool> vibSpecieActive(Y_species.size(), false);

// No change from original — vibSpecieActive[N] = false
// vibEnergyEquation.H skips FVM solve for atoms
// updateVibTemperature.H handles atomic electronic diagnostically
forAll(Y_species, speciei)
{
    vibSpecieActive[speciei] =
        solveSpeciesVib && mix.isSpecieMolecular(speciei);
}

// sensible enthalpy
volScalarField& hTR = thermo2T.h();
volScalarField& eTR = thermo2T.eTR();
volScalarField& eV = thermo2T.eVib();

if (use2T)
{
    eV.correctBoundaryConditions();
}
else
{
    eV = dimensionedScalar("zero", dimEnergy/dimMass, 0.0);

    eV.correctBoundaryConditions();
}

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
    thermo2T.rho()
);

PtrList<volScalarField> TVibSpecies;
PtrList<volScalarField> eVibSpecies;
PtrList<volScalarField> rhoYEvibSpecies;

autoPtr<volScalarField> rhoEvPtr
(
    new volScalarField
    (
        IOobject
        (
            "rhoEv",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume, 0.0)
    )
);
volScalarField& rhoEv = rhoEvPtr();

if (use2T)
{
    rhoEv = rho*eV;
}
else
{
    rhoEv = dimensionedScalar("zero", dimEnergy/dimVolume, 0.0);
}
rhoEv.correctBoundaryConditions();

if (solveSpeciesVib)
{
    TVibSpecies.setSize(Y_species.size());
    eVibSpecies.setSize(Y_species.size());
    rhoYEvibSpecies.setSize(Y_species.size());

    forAll(Y_species, speciei)
    {
        const word sName = Y_species[speciei].name();

        TVibSpecies.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "TVib_" + sName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                TVib
            )
        );

        eVibSpecies.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "eVib_" + sName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                eV
            )
        );

        rhoYEvibSpecies.set
        (
            speciei,
            new volScalarField
            (
                IOobject
                (
                    "rhoYEvib_" + sName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimEnergy/dimVolume, 0.0)
            )
        );
    }

    const scalarField& pCells       = p.primitiveField();
    const scalarField& TTRCells     = TTR.primitiveField();
    const scalarField& TVibMixCells = TVib.primitiveField();

    const volScalarField::Boundary& pBf       = p.boundaryField();
    const volScalarField::Boundary& TTRBf     = TTR.boundaryField();
    const volScalarField::Boundary& TVibMixBf = TVib.boundaryField();
    
    const scalarField& rhoCells = rho.primitiveField();
    const volScalarField::Boundary& rhoBf = rho.boundaryField();

    forAll(Y_species, speciei)
    {
        volScalarField& TVs  = TVibSpecies[speciei];
        volScalarField& eVs  = eVibSpecies[speciei];
        volScalarField& rYeV = rhoYEvibSpecies[speciei];

        if (!vibSpecieActive[speciei])
        {
            TVs  = dimensionedScalar("zero", dimTemperature, 0.0);
            eVs  = dimensionedScalar("zero", dimEnergy/dimMass, 0.0);
            rYeV = dimensionedScalar("zero", dimEnergy/dimVolume, 0.0);

            TVs.correctBoundaryConditions();
            eVs.correctBoundaryConditions();
            rYeV.correctBoundaryConditions();
            continue;
        }

        const scalar thetaVib_i = mix.thetaVib(speciei);

        TVs = TVib;

        scalarField& eVsCells  = eVs.primitiveFieldRef();
        scalarField& rYeVCells = rYeV.primitiveFieldRef();
        const scalarField& YiCells = Y_species[speciei].primitiveField();

        forAll(eVsCells, celli)
        {
            eVsCells[celli] =
                mix.EsVib
                (
                    speciei,
                    pCells[celli],
                    TTRCells[celli],
                    TVibMixCells[celli],
                    thetaVib_i
                );

            rYeVCells[celli] = rhoCells[celli]*YiCells[celli]*eVsCells[celli];
        }

        volScalarField::Boundary& eVsBf  = eVs.boundaryFieldRef();
        volScalarField::Boundary& rYeVBf = rYeV.boundaryFieldRef();
        const volScalarField::Boundary& YiBf = Y_species[speciei].boundaryField();

        forAll(pBf, patchi)
        {
            const fvPatchScalarField& pp    = pBf[patchi];
            const fvPatchScalarField& pTTR  = TTRBf[patchi];
            const fvPatchScalarField& pTVib = TVibMixBf[patchi];
            const fvPatchScalarField& prho  = rhoBf[patchi];
            const fvPatchScalarField& pYi   = YiBf[patchi];

            fvPatchScalarField& peVib = eVsBf[patchi];
            fvPatchScalarField& prYeV = rYeVBf[patchi];

            forAll(pp, facei)
            {
                peVib[facei] =
                    mix.EsVib
                    (
                        speciei,
                        pp[facei],
                        pTTR[facei],
                        pTVib[facei],
                        thetaVib_i
                    );

                prYeV[facei] = prho[facei]*pYi[facei]*peVib[facei];
            }
        }

        TVs.correctBoundaryConditions();
        eVs.correctBoundaryConditions();
        rYeV.correctBoundaryConditions();
    }
}

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

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(eTR + eV + 0.5*magSqr(U))
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info<< "Turbulence model: Laminar\n" << endl;


// Mixture rho for sanity check purposes
volScalarField rhoCont
(
    IOobject("rhoCont", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
    rho
);

// Mixture fluxes
PtrList<volScalarField> rhoYi;

if (mixtureCheck)
{
    PtrList<volScalarField>& Y = thermo2T.composition().Y();
    rhoYi.setSize(Y.size());

    forAll(Y, i)
    {
        const word name("rhoY_" + Y[i].name());

        rhoYi.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    name,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rho*Y[i]
            )
        );
    }
}

// Effective diffusion coefficient field [m^2/s]
volScalarField Deff
(
    IOobject
    (
        "Deff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE     
    ),
    mesh,
    dimensionedScalar
    (
        "zero", 
        dimArea/dimTime, 
        0.0
    )
);

// --- Fourier conductive heat flux for Eq. (7): Translational-Rotational Energy Equation
surfaceScalarField qFourierPhi
(
    IOobject
    (
        "qFourierPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimEnergy/dimTime,
        0.0
    )
);

PtrList<surfaceScalarField> phiJYi;

if (mixtureCheck)
{
    const PtrList<volScalarField>& Y = thermo2T.composition().Y();

    phiJYi.setSize(Y.size());

    forAll(Y, i)
    {
        phiJYi.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phiJYi_" + Y[i].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimMass/dimTime,
                    0.0
                )
            )
        );

        phiJYi[i].setOriented(true);
    }
}

surfaceScalarField qVibPhi
(
    IOobject
    (
        "qVibPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimEnergy/dimTime,
        0.0
    )
);

qVibPhi.setOriented(true);
