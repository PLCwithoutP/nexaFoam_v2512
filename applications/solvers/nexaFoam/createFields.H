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

// ---- P3-06 (A2): 2T needs a multi-component mixture, also for one species
//      (as D7). The pure-mixture 2T path (thermo correctTVib, A2) is not used.
if (use2T && !mixtureCheck)
{
    FatalErrorInFunction
        << "twoTemperature true requires a multi-component mixture "
        << "(multiComponent2TMixture or neReactingMixture), also for a "
        << "single species." << nl
        << exit(FatalError);
}

// ---- P3-01 (G1): start-up check that the closed-form inverse TvHO is the
//      exact inverse of the library EsVib for every molecule.
if (solveSpeciesVib)
{
    const scalar Ttest[4] = {50, 300, 3000, 20000};

    forAll(Y_species, speciei)
    {
        if (!vibSpecieActive[speciei]) continue;

        const scalar theta = mix.thetaVib(speciei);
        const scalar R     = mix.R(speciei);

        for (label k = 0; k < 4; ++k)
        {
            const scalar ev =
                mix.EsVib(speciei, 1e5, Ttest[k], Ttest[k], theta);
            const scalar Tback = vibHelpers::TvHO(ev, R, theta);

            if (!(mag(Tback - Ttest[k]) <= 1e-10*Ttest[k]))
            {
                FatalErrorInFunction
                    << "G1 check failed for " << Y_species[speciei].name()
                    << ": T = " << Ttest[k] << " K, EsVib = " << ev
                    << " J/kg, closed-form inverse = " << Tback << " K." << nl
                    << "The library vibrational model is no longer the "
                    << "harmonic oscillator with R = mix.R(speciei)."
                    << exit(FatalError);
            }
        }
    }

    Info<< "G1 check passed: the closed-form T_v inverse matches EsVib "
        << "for every molecule (50, 300, 3000, 20000 K)." << nl << endl;
}

// C-V model data for the chemistry split (applyChemistry.H)
// The same keywords the library CVEnergySource reads, so the cases need no
// change: CVTreatment (preferential | non-preferential, required) and
// alphaDissociation (optional, default 0.3) in thermophysicalProperties,
// dissociationPotential [J/kg] per molecule in speciesDict.
word CVTreatment("non-preferential");
scalar alphaDissociation = 0.3;
scalarList Dissoc(Y_species.size(), Zero);

if (use2T && mixtureCheck && useChemistry)
{
    const IOdictionary thermoDictCV
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    CVTreatment = thermoDictCV.get<word>("CVTreatment");

    if (CVTreatment != "preferential" && CVTreatment != "non-preferential")
    {
        FatalIOErrorInFunction(thermoDictCV)
            << "Invalid CVTreatment = '" << CVTreatment << "'." << nl
            << "Valid values are 'preferential' and 'non-preferential'."
            << exit(FatalIOError);
    }

    thermoDictCV.readIfPresent("alphaDissociation", alphaDissociation);

    const IOdictionary speciesDictCV
    (
        IOobject
        (
            "speciesDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    forAll(Y_species, speciei)
    {
        if (!vibSpecieActive[speciei]) continue;

        // Field name "Y_<specie>" -> speciesDict key "<specie>"
        const word fName(Y_species[speciei].name());
        const word sName
        (
            fName.substr(0, 2) == "Y_" ? word(fName.substr(2)) : fName
        );

        Dissoc[speciei] =
            speciesDictCV.subDict(sName).get<scalar>("dissociationPotential");
    }

    Info<< "C-V exchange in the chemistry split: " << CVTreatment;
    if (CVTreatment == "preferential")
    {
        Info<< ", alphaDissociation = " << alphaDissociation;
    }
    Info<< nl << endl;
}

const bool CVpreferential = (CVTreatment == "preferential");
// ---- end P5-01 -----------------------------------------------------------------

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

{
    const label cWall = mesh.boundaryMesh().findPatchID("cone");   // your wall patch name
    if (cWall >= 0)
    {
        Info<< "[BCcheck] cone eV   type = " << eV.boundaryField()[cWall].type()   << nl
            << "[BCcheck] cone TVib type = " << TVib.boundaryField()[cWall].type() << endl;
    }
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

        // A per-molecule T_v field is read from 0/TVib_<Y name>
        // when the case provides one (per-molecule wall BC, e.g.
        // neSmoluchowskiJumpTv); otherwise it is a copy of the mixture TVib.
        IOobject TVsIO
        (
            "TVib_" + sName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (TVsIO.typeHeaderOk<volScalarField>(true))
        {
            Info<< "Reading per-molecule field " << TVsIO.name() << nl;
            TVibSpecies.set(speciei, new volScalarField(TVsIO, mesh));
        }
        else
        {
            TVsIO.readOpt(IOobject::NO_READ);
            TVibSpecies.set(speciei, new volScalarField(TVsIO, TVib));
        }

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

    eV = dimensionedScalar("zero", eV.dimensions(), 0.0);
    forAll(Y_species, speciei)
    {
        if (!vibSpecieActive[speciei]) continue;
        eV += rhoYEvibSpecies[speciei]/rho;
    }
    eV.correctBoundaryConditions();

    rhoEv = rho*eV;
    rhoEv.correctBoundaryConditions();
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

volScalarField muLam
(
    IOobject
    (
        "muLam",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo2T.mu()
);

volScalarField kappaPP
(
    IOobject
    (
        "kappaPP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo2T.kappaTR()          
);

volScalarField kappaVe
(
    IOobject
    (
        "kappaVe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.0*kappaPP            // inherits W/(m·K) dimensions, initialised to zero
);

if (use2T)
{
    kappaVe = thermo2T.kappaVib();   
}

// ---- Transport, per-molecule jump and
//      wall heat-flux fields ---------------------------------------------------

// Lewis number of the Fick model (C5): keyword Le in thermophysicalProperties,
// default 1.4 as before
const scalar LeNum =
    IOdictionary
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    ).getOrDefault<scalar>("Le", 1.4);

Info<< "Fick diffusion: Le = " << LeNum << nl << endl;

// Wilke-mixed kappa_TR (transportUpdate.H); registered: the TR jump BC
// (neSmoluchowskiJumpTR) reads its patch values
volScalarField kappaTRW
(
    IOobject
    (
        "kappaTRW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo2T.kappaTR()
);

// Frozen TR heat capacity of the mixture, c_v,TR = sum_s Y_s (c_p,s - R_s)
volScalarField cvTRmix
(
    IOobject
    (
        "cvTRmix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    ),
    mesh,
    dimensionedScalar("one", dimEnergy/dimMass/dimTemperature, 1.0)
);

// Gamma_v,m = Y_m Gh_m (G3), read by the per-molecule T_v jump BC
PtrList<volScalarField> GammaVibSpecies(Y_species.size());
forAll(Y_species, speciei)
{
    GammaVibSpecies.set
    (
        speciei,
        new volScalarField
        (
            IOobject
            (
                "GammaVib_" + Y_species[speciei].name(),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimMass/dimLength/dimTime, 0.0)
        )
    );
}

// Wall heat flux from the solver's own conduction fluxes [W/m^2],
// positive into the wall (D4.5)
volScalarField qWallTR
(
    IOobject
    (
        "qWallTR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0)
);

volScalarField qWallVib
(
    IOobject
    (
        "qWallVib",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimTime/dimArea, 0.0)
);

kappaVe.correctBoundaryConditions();

// --- Positivity / admissibility controls (system/fvSolution: positivity{})
const scalar positivityTFloor =
    mesh.solutionDict().subOrEmptyDict("positivity")
        .getOrDefault<scalar>("TFloor", 25.0);

Info<< "Positivity: TFloor = " << positivityTFloor << " K" << nl << endl;