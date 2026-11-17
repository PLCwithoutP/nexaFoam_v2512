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

{
    const scalar TFloor = positivityTFloor;
    const scalarField& V        = mesh.V();
    const scalarField& rhoCells = rho.primitiveField();
    const scalarField& eVCells  = eV.primitiveField();
    const scalarField& pCells   = p.primitiveField();

    scalarField& rhoECells = rhoE.primitiveFieldRef();
    vectorField& rhoUCells = rhoU.primitiveFieldRef();
    vectorField& UCells    = U.primitiveFieldRef();

    label  nRescaled = 0;
    label  nReset    = 0;
    scalar dMomentum = 0;   // |momentum| removed  [kg m/s]
    scalar dEnergy   = 0;   // energy added        [J]

    forAll(rhoCells, celli)
    {
        const scalar rhoc = rhoCells[celli];
        if (!(rhoc > 0)) continue;          // caught by the density check

        // e_TR of this cell's composition at TFloor
        scalar eMin = 0;
        forAll(Y_species, speciei)
        {
            eMin +=
                Y_species[speciei][celli]
               *(
                    mix.Hs(speciei, pCells[celli], TFloor)
                  - mix.R(speciei)*TFloor
                );
        }

        const scalar eAvail = rhoECells[celli]/rhoc - eVCells[celli] - eMin;
        const scalar KE     = 0.5*magSqr(UCells[celli]);

        if (KE <= eAvail) continue;         // admissible: nothing to do

        if (eAvail > 0)
        {
            // Scale velocity so that KE = eAvail -> e_TR = eMin exactly
            const scalar s = Foam::sqrt(eAvail/KE);
            dMomentum += rhoc*V[celli]*mag(UCells[celli])*(1.0 - s);
            UCells[celli]   *= s;
            rhoUCells[celli] = rhoc*UCells[celli];
            ++nRescaled;
        }
        else
        {
            // Not enough energy even at rest: stop the cell, lift e_TR
            dMomentum += rhoc*V[celli]*mag(UCells[celli]);
            const scalar rhoENew = rhoc*(eMin + eVCells[celli]);
            dEnergy += (rhoENew - rhoECells[celli])*V[celli];
            UCells[celli]    = Zero;
            rhoUCells[celli] = Zero;
            rhoECells[celli] = rhoENew;
            ++nReset;
        }
    }

    reduce(nRescaled, sumOp<label>());
    reduce(nReset,    sumOp<label>());
    reduce(dMomentum, sumOp<scalar>());
    reduce(dEnergy,   sumOp<scalar>());

    if (nRescaled + nReset > 0)
    {
        Info<< "    [admissibility] e_TR < e_TR(" << TFloor << " K): "
            << nRescaled << " cells KE-rescaled, "
            << nReset    << " cells reset; momentum removed = " << dMomentum
            << " kg m/s, energy added = " << dEnergy << " J" << nl;
    }

    U.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();
}

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

    // TR conduction with the Wilke kappa_TR,
    // implicit on e_TR with Gamma_TR = kappa_TR,W/c_v,TR and the composition
    // term of the chain rule,
    //     kappa grad T = Gamma_TR (grad e_TR - sum_s e_TR,s grad Y_s),
    // exact also where c_v,TR varies (molecules 5/2 R, atoms 3/2 R).
    // Pure mixtures keep the library form (alphaheTR*gamma).
    // Explicit: species diffusive enthalpy, viscous work.
    volScalarField GammaTR
    (
        "GammaTR",
        mixtureCheck
    ? kappaTRW/cvTRmix
    : thermo2T.alphaheTR()*thermo2T.gamma()
    );

    surfaceScalarField phiTRcomp
    (
        IOobject("phiTRcomp", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    );
    phiTRcomp.setOriented(true);

    if (mixtureCheck)
    {
        forAll(Y_species, speciei)
        {
            const scalar Rs = mix.R(speciei);

            volScalarField eTRs
            (
                IOobject("eTR_" + Y_species[speciei].name(), runTime.timeName(),
                        mesh, IOobject::NO_READ, IOobject::NO_WRITE,
                        IOobject::NO_REGISTER),
                mesh,
                dimensionedScalar("zero", eTR.dimensions(), 0.0)
            );

            forAll(eTRs, celli)
            {
                eTRs[celli] =
                    mix.Hs(speciei, p[celli], TTR[celli]) - Rs*TTR[celli];
            }

            forAll(eTRs.boundaryField(), patchi)
            {
                fvPatchScalarField& pe = eTRs.boundaryFieldRef()[patchi];
                const fvPatchScalarField& pp = p.boundaryField()[patchi];
                const fvPatchScalarField& pT = TTR.boundaryField()[patchi];

                forAll(pe, facei)
                {
                    pe[facei] = mix.Hs(speciei, pp[facei], pT[facei]) - Rs*pT[facei];
                }
            }

            phiTRcomp +=
                fvc::interpolate(GammaTR*eTRs)
                *fvc::snGrad(Y_species[speciei])*mesh.magSf();
        }
    }

    solve
    (
        fvm::ddt(rho, eTR) - fvc::ddt(rho, eTR)
    - fvm::laplacian(GammaTR, eTR)
    ==
    - fvc::div(phiJhTR + phiJhVib)
    - fvc::div(phiTRcomp)
    );

    // ---- P4-07 (D4.5): TR wall heat flux from the discrete flux of this
    // solve, positive into the wall
    {
        const surfaceScalarField GammaTRf(fvc::interpolate(GammaTR));
        const surfaceScalarField sgE(fvc::snGrad(eTR));

        forAll(mesh.boundary(), patchi)
        {
            if (!isA<wallFvPatch>(mesh.boundary()[patchi])) continue;

            qWallTR.boundaryFieldRef()[patchi] =
                phiTRcomp.boundaryField()[patchi]
                /mesh.magSf().boundaryField()[patchi]
            - GammaTRf.boundaryField()[patchi]*sgE.boundaryField()[patchi];
        }

        if (debugLevel >= 1)
        {
            scalar QTR = 0;
            forAll(mesh.boundary(), patchi)
            {
                if (!isA<wallFvPatch>(mesh.boundary()[patchi])) continue;
                QTR += sum
                (
                    qWallTR.boundaryField()[patchi]
                    *mesh.magSf().boundaryField()[patchi]
                );
            }
            Info<< "  [qWall] TR = " << returnReduce(QTR, sumOp<scalar>())
                << " W" << nl;
        }
    }
    // ---- end P4-04 / P4-07 -----------------------------------------------------

    thermo2T.correct();
    // Sync rhoE with viscously corrected eTR
    rhoE = rho*(eTR + eV + 0.5*magSqr(U));    
}

p.ref() = rho() / psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
