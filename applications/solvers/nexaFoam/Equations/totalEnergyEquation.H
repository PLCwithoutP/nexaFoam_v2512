// ================================================================
// Stage 0 — Viscous Work
// ================================================================
surfaceScalarField sigmaDotU
(
    "sigmaDotU",
    (
        fvc::interpolate(muLam)*mesh.magSf()*fvc::snGrad(U)
        + fvc::dotInterpolate(mesh.Sf(), tauMC)
    )
    & (kurganovFlux.alpha()*kurganovFlux.U_pos() + kurganovFlux.minus_alpha()*kurganovFlux.U_neg())
);

// ================================================================
// Stage 1 — Inviscid explicit
// ================================================================
solve
(
    fvm::ddt(rhoE)
  + fvc::div(phiEp)
  - fvc::div(sigmaDotU)
);

// ================================================================
// Stage 2 — Primitive update  [absorbs updateTTRTemperature.H]
// ================================================================
thermo2T.eVib() = rhoEv / rho;

eTR =
    rhoE / rho
  - 0.5*magSqr(U)
  - eV;

eTR.correctBoundaryConditions();
thermo2T.correct();

rhoE.boundaryFieldRef() ==
    rho.boundaryField()
   *(
        eTR.boundaryField()
      + eV.boundaryField()
      + 0.5*magSqr(U.boundaryField())
    );



// ================================================================
// Stage 3 — Viscous correction
// ================================================================
if (!inviscid)
{
    // Species diffusive enthalpy flux (explicit)
    surfaceScalarField phiJhTR
    (
        IOobject("phiJhTR", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    );
    phiJhTR.setOriented(true);

    surfaceScalarField phiJhVib
    (
        IOobject("phiJhVib", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    );
    phiJhVib.setOriented(true);

    if (mixtureCheck)
    {
        forAll(Y_species, speciei)
        {
            volScalarField hsTRi
            (
                IOobject("hsTR_" + Y_species[speciei].name(), runTime.timeName(),
                         mesh, IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
                mesh,
                dimensionedScalar("zero", hTR.dimensions(), 0.0)
            );

            forAll(hsTRi, celli)
            {
                hsTRi[celli] = mix.Hs(speciei, p[celli], TTR[celli]);
            }

            forAll(hsTRi.boundaryField(), patchi)
            {
                forAll(hsTRi.boundaryField()[patchi], facei)
                {
                    hsTRi.boundaryFieldRef()[patchi][facei] =
                        mix.Hs(speciei,
                               p.boundaryField()[patchi][facei],
                               TTR.boundaryField()[patchi][facei]);
                }
            }

            phiJhTR += phiJYi[speciei]*fvc::interpolate(hsTRi);
        }
    }

    if (mixtureCheck && use2T)
    {
        forAll(Y_species, speciei)
        {
            if (!vibSpecieActive[speciei]) continue;
            phiJhVib += phiJYi[speciei]*fvc::interpolate(eVibSpecies[speciei]);
        }
    }

    // Implicit TR conduction via alphaheTR * gamma = kappaTR/CvTR
    // Explicit: species diffusive enthalpy, viscous work
    volScalarField alphavTR
    (
        "alphavTR",
        thermo2T.alphaheTR() * thermo2T.gamma()
    );

    solve
    (
        fvm::ddt(rho, eTR) - fvc::ddt(rho, eTR)
        - fvm::laplacian(alphavTR, eTR)
     ==
        - fvc::div(phiJhTR + phiJhVib)
    );

    thermo2T.correct();
    // Sync rhoE with viscously corrected eTR
    rhoE = rho*(eTR + eV + 0.5*magSqr(U));    
}

p.ref() = rho() / psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
