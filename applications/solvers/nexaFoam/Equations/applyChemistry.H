if (useChemistry)
{
    const scalar dtChem = runTime.deltaTValue();

    // 1) Solve the local chemistry ODE system.
    //    This internally uses the chemistry solver selected in chemistryProperties (currently neEulerImplicit).
    chemistryPtr->solve(dtChem);

    // 2) Species source-term split update: rhoYi^{n+1} = rhoYi^* + dt * RR_i
    PtrList<volScalarField>& Y = thermo2T.composition().Y();
    const label nSpec = Y.size();

    label nNegChem = 0;     // P3-05 (D3.2)
    for (label i = 0; i < nSpec; ++i)
    {
        const volScalarField::Internal& RRi = chemistryPtr->RR(i);

        forAll(rhoYi[i], celli)
        {
            rhoYi[i][celli] += dtChem*RRi[celli];

            // P3-05 (D3.2): no floor on the conserved rhoYi. An overshoot of
            // the split chemistry update is counted (B4, Phase 6) and clipped
            // only in the property field Y (step 4 below).
            if (rhoYi[i][celli] < -1e-12*rho[celli])
            {
                ++nNegChem;
            }
        }

        rhoYi[i].correctBoundaryConditions();
    }
    reduce(nNegChem, sumOp<label>());
    if (nNegChem > 0)
    {
        Info<< "    [chem] rhoY_s < -1e-12 rho after the chemistry update kept: "
            << nNegChem << " species-cells" << nl;
    }

    // 3) Rebuild mixture density from partial densities
    rho = dimensionedScalar("zeroRho", rho.dimensions(), 0.0);

    for (label i = 0; i < nSpec; ++i)
    {
        rho += rhoYi[i];
    }

    rho.max(SMALL);
    rho.correctBoundaryConditions();

    // 4) Recover mass fractions Yi = rhoYi / rho
    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] = rhoYi[i]/rho;
        Y[i].max(0.0);
        Y[i].min(1.0);
        Y[i].correctBoundaryConditions();
    }

    // 5) Renormalize Y
    volScalarField sumY
    (
        IOobject
        (
            "sumYChemistry",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    );

    for (label i = 0; i < nSpec; ++i)
    {
        sumY += Y[i];
    }

    sumY.max(SMALL);

    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] /= sumY;
        Y[i].correctBoundaryConditions();
    }

    // Chemistry-vibration exchange in the same split step as the species
    // update and with the same d(rhoY_m) = dt*omegaDot_m (no lag, K8).
    // Destruction (d(rhoY_m) < 0):
    //   non-preferential: the molecules removed carry the molecule's own
    //       e_v,m = eps_m/rhoY_m (I3), i.e. e_v,m is unchanged:
    //       eps_new = eps*rhoY_new/rhoY_old  (exact, >= 0 while rhoY_new >= 0)
    //   preferential: they carry alphaDissociation*D_m, Patankar-implicit:
    //       eps_new = eps^2/(eps + |d rhoY| alpha D_m)  (>= 0 for any dt, G6)
    // Formation (d(rhoY_m) > 0, D5.4): the new molecules carry e_v,m(T_v,m)
    //   of the cell (non-preferential) or alphaDissociation*D_m (preferential).
    // The energy only moves between eps_m and e_TR: rhoE is not touched, and
    // e_TR is rebuilt from rhoE in step 6/7 with the eV computed here.
    if (use2T && mixtureCheck)
    {
        label  nPatankar = 0;   // cells where alpha*D_m*|d rhoY| > eps
        scalar dEcv      = 0;   // energy moved into vibration [J]
        const scalarField& Vc   = mesh.V();
        const scalarField& rhoC = rho.primitiveField();
        const scalarField& TTRC = TTR.primitiveField();

        forAll(Y_species, m)
        {
            if (!vibSpecieActive[m]) continue;

            const scalarField& RRm = chemistryPtr->RR(m);
            const scalarField& rYm = rhoYi[m].primitiveField();
            const scalarField& evm = eVibSpecies[m].primitiveField();
            scalarField& eps       = rhoYEvibSpecies[m].primitiveFieldRef();
            const scalar eHatPref  = alphaDissociation*Dissoc[m];

            forAll(eps, celli)
            {
                const scalar dRhoY = dtChem*RRm[celli];
                if (dRhoY == 0) continue;

                const scalar eps0    = eps[celli];
                const scalar rhoYnew = rYm[celli];          // after step 2
                const scalar rhoYold = rhoYnew - dRhoY;

                if (dRhoY > 0)
                {
                    eps[celli] +=
                        dRhoY
                        *(CVpreferential ? eHatPref : max(evm[celli], scalar(0)));
                }
                else if (CVpreferential)
                {
                    if (eps0 > VSMALL)
                    {
                        const scalar loss = -dRhoY*eHatPref;
                        if (loss > eps0) ++nPatankar;
                        eps[celli] = eps0*eps0/(eps0 + loss);
                    }
                }
                else if (rhoYold > VSMALL)
                {
                    eps[celli] = eps0*max(rhoYnew, scalar(0))/rhoYold;
                }

                dEcv += (eps[celli] - eps0)*Vc[celli];
            }

            rhoYEvibSpecies[m].correctBoundaryConditions();

            eVibSpecies[m] =
                rhoYEvibSpecies[m]/max(rhoYi[m], vibHelpers::YTrace*rho);
            eVibSpecies[m].correctBoundaryConditions();

            // T_v,m of the molecules present (G1, D3.3); trace molecules keep
            // theirs until the next Stage 2 (G4)
            scalarField& TVs        = TVibSpecies[m].primitiveFieldRef();
            const scalarField& eVs  = eVibSpecies[m].primitiveField();
            const scalar Rm         = mix.R(m);
            const scalar thm        = mix.thetaVib(m);
            forAll(TVs, celli)
            {
                if (rYm[celli] < vibHelpers::YTrace*rhoC[celli]) continue;
                TVs[celli] =
                    eVs[celli] > VSMALL
                ? vibHelpers::TvHO(eVs[celli], Rm, thm)
                : TTRC[celli];
            }
            TVibSpecies[m].correctBoundaryConditions();
        }

        // Mixture e_v from the new eps and rho, used for e_TR in step 6/7
        eV = dimensionedScalar("zero", eV.dimensions(), 0.0);
        forAll(Y_species, m)
        {
            if (!vibSpecieActive[m]) continue;
            eV += rhoYEvibSpecies[m]/rho;
        }
        eV.correctBoundaryConditions();

        Info<< "    [CV] " << CVTreatment << ": energy into vibration = "
            << returnReduce(dEcv, sumOp<scalar>()) << " J, Patankar-limited = "
            << returnReduce(nPatankar, sumOp<label>()) << " molecule-cells" << nl;
    }

    // Step 6 & 7: Energy update from chemistry source
    //
    // neTCLib:    hTR = sensible enthalpy only → Y-change gives ΔHs,
    //             Qdot = -sum(Hf_i*RR_i) provides the missing ΔHf.  NEEDED.
    //
    // Mutation++: hTR = mixtureHMass (total, includes Hf) → Y-change gives
    //             ΔHs+ΔHf = full physical cooling already.  Applying Qdot
    //             here double-counts ΔHf and causes ~2x excess cooling.  SKIP.
    if (thermo2T.chemistryReader() != word("mutationpp"))
    {
        // neTCLib path: hTR = sensible enthalpy only.
        // Y-change provides ΔHs; Qdot provides the missing ΔHf.
        tmp<volScalarField> tQdot = chemistryPtr->Qdot();
        const volScalarField& Qdot = tQdot();

        forAll(rhoE, celli)
        {
            rhoE[celli] += dtChem * Qdot[celli];
        }
        rhoE.correctBoundaryConditions();

        eTR =
            rhoE/rho
        - 0.5*magSqr(U)
        - eV;
        eTR.correctBoundaryConditions();
    }
    else
    {
        // Mutation++ path: hTR = mixtureHMass (total enthalpy, includes Hf).
        // The Y-change from steps 1–5 already encodes the full thermodynamic
        // cooling: the next updateTTRTemperature will invert rhoE against the
        // new composition and find a physically correct lower TTR.
        // Applying Qdot to rhoE here would double-count the formation enthalpy
        // contribution and cause ~2x excess cooling.
        // No explicit rhoE or hTR modification is required.
    }

    // Keep momentum consistent with unchanged velocity
    rhoU = rho*U;
    rhoU.correctBoundaryConditions();

    // 8) Re-correct thermo from the chemistry-updated state
    thermo2T.correct();

    if (use2T && mixtureCheck)
    {
        // Rebuild eV from per-species fields — consistent with vibEnergyEquation
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
    else if (use2T)
    {
        thermo2T.correctVibEnergy();
    }

    // 9) Rebuild pressure from updated rho and psi
    p.ref() = rho()/psi();
    p.correctBoundaryConditions();
    rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
    
    // Optional consistency for mixture-check path
    if (mixtureCheck)
    {
        rhoCont = rho;
    }
}
