if (use2T)
{
    volScalarField QVzero
    (
        IOobject("QVzero", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );

    volScalarField pe_divU
    (
        IOobject("pe_divU", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        thermo2T.pe()() * fvc::div(U)
    );

    // ============================================================
    // Stage 1 — Inviscid + sources
    // ============================================================
    if (!mixtureCheck)
    {
        solve
        (
            fvm::ddt(rhoEv)
          + fvc::div(phiEv)
         ==
            QVzero
        );

        eV = rhoEv/rho;
        eV.correctBoundaryConditions();
    }
    else
    {
        PtrList<volScalarField>& QVT = thermo2T.correctVibSource(TVibSpecies);

        PtrList<volScalarField> QVVempty;
        PtrList<volScalarField>* QVVptr = &QVVempty;
        if (useVibVibRelaxation)
        {
            QVVptr = &thermo2T.correctVibVibSource(TVibSpecies);
        }
        PtrList<volScalarField>& QVV = *QVVptr;

        // The C-V exchange is applied in applyChemistry.H, in
        // the same split step as the species update.

        PtrList<volScalarField>& Y = thermo2T.composition().Y();

        forAll(Y, speciei)
        {
            if (!vibSpecieActive[speciei])
            {
                rhoYEvibSpecies[speciei] =
                    dimensionedScalar("zero", rhoYEvibSpecies[speciei].dimensions(), 0.0);
                eVibSpecies[speciei] =
                    dimensionedScalar("zero", eVibSpecies[speciei].dimensions(), 0.0);
                TVibSpecies[speciei] =
                    dimensionedScalar("zero", TVibSpecies[speciei].dimensions(), 0.0);
                rhoYEvibSpecies[speciei].correctBoundaryConditions();
                eVibSpecies[speciei].correctBoundaryConditions();
                TVibSpecies[speciei].correctBoundaryConditions();
                continue;
            }

            surfaceScalarField phiJEvib_i
            (
                IOobject("phiJEvib_" + Y[speciei].name(), runTime.timeName(),
                         mesh, IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
                phiJYi[speciei]*fvc::interpolate(eVibSpecies[speciei])
            );
            phiJEvib_i.setOriented(true);

            const volScalarField& QVTi = (speciei < QVT.size() ? QVT[speciei] : QVzero);
            const volScalarField& QVVi = (speciei < QVV.size() ? QVV[speciei] : QVzero);

            solve
            (
                fvm::ddt(rhoYEvibSpecies[speciei])
              + fvc::div(phiRhoEvibYi[speciei] + phiJEvib_i)
             ==
                QVTi + QVVi
            );

            // ---- P3-02 (G4, D3.2) ------------------------------------------
            // No floor on the conserved rhoYi, no clip of eps_m = rhoYEvib.
            // e_v,m = eps_m/rhoY_m with the denominator bounded below by
            // YTrace*rho, so a trace molecule never gets an absurd e_v (K2);
            // its cells are reset by the trace guard in Stage 2. A negative
            // eps_m is kept (conservation) and gives T_v = T_TR (D3.3).
            rhoYEvibSpecies[speciei].correctBoundaryConditions();

            eVibSpecies[speciei] =
                rhoYEvibSpecies[speciei]
                /max(rhoYi[speciei], vibHelpers::YTrace*rho);
            eVibSpecies[speciei].correctBoundaryConditions();
            // ---- end P3-02 -------------------------------------------------
        }

        //  eps_m < 0 after Stage 1 (convection,
        //  Fick, V-T, V-V), before any decode. Counted, not clipped.
        {
            label nNegEps = 0;
            forAll(Y, speciei)
            {
                if (!vibSpecieActive[speciei]) continue;
                for (const scalar e : rhoYEvibSpecies[speciei].primitiveField())
                {
                    if (e < 0) ++nNegEps;
                }
            }
            reduce(nNegEps, sumOp<label>());
            if (nNegEps > 0)
            {
                Info<< "    [vib] eps < 0 after Stage 1: " << nNegEps
                    << " molecule-cells" << nl;
            }
        }

        eV = dimensionedScalar("zero", eV.dimensions(), 0.0);
        forAll(Y, speciei)
        {
            if (!vibSpecieActive[speciei]) continue;
            eV += rhoYEvibSpecies[speciei]/rho;
        }
        eV.correctBoundaryConditions();
        rhoEv = rho*eV;
        rhoEv.correctBoundaryConditions();
        rhoEv -= runTime.deltaT() * pe_divU;
        rhoEv.correctBoundaryConditions();
        eV = rhoEv/rho;
        eV.correctBoundaryConditions();
    }

    // ============================================================
    // Stage 2 — Primitive update
    // ============================================================
    if (!mixtureCheck)
    {
        thermo2T.correctTVib();
        TVib.correctBoundaryConditions();
        eV.correctBoundaryConditions();
        rhoEv = rho*eV;
        rhoEv.correctBoundaryConditions();
    }
    else
    {
        // ---- P3-03 (G1, G4, G2, D3.3) ---------------------------------------
        // Per molecule: exact closed-form inversion (G1) where Y_m >= YTrace;
        // T_v = T_TR where e_v <= 0 (D3.3); trace molecules take the mixture T_v
        // (G4). Mixture T_v: energy-equivalent (G2). Replaces the e_v floor of
        // 20 J/kg (K1), the library Newton (A6) and the mass-weighted mean (K4).
        // No conserved field (rhoYi, rhoYEvibSpecies) is modified here.
        const scalarField& pCells   = p.primitiveField();
        const scalarField& TTRCells = TTR.primitiveField();
        const scalarField& rhoCells = rho.primitiveField();
        scalarField& TVibMixCells   = TVib.primitiveFieldRef();

        const label nSp = Y_species.size();
        scalarList RSp(nSp, Zero);
        scalarList thetaSp(nSp, Zero);
        forAll(Y_species, speciei)
        {
            if (!vibSpecieActive[speciei]) continue;
            RSp[speciei]     = mix.R(speciei);
            thetaSp[speciei] = mix.thetaVib(speciei);
        }

        // Work lists for the G2 Newton (molecules present in one cell or face)
        scalarList wY(nSp), wR(nSp), wTh(nSp), wTv(nSp), wEv(nSp);

        label  nTrace   = 0;    // molecule-cells set by the trace guard (G4)
        label  nNonPos  = 0;    // molecule-cells with e_v <= 0: T_v = T_TR (D3.3)
        scalar maxRT    = 0;    // G1 round trip: max |EsVib(T_v) - e_v|/e_v
        scalar maxG2    = 0;    // G2: max relative energy mismatch
        scalar maxG2f   = 0;    // same on boundary faces (not reported)

        // (a) Cells, molecules present: G1, or T_TR for e_v <= 0
        forAll(Y_species, speciei)
        {
            if (!vibSpecieActive[speciei])
            {
                TVibSpecies[speciei] =
                    dimensionedScalar("zero", dimTemperature, 0.0);
                TVibSpecies[speciei].correctBoundaryConditions();
                continue;
            }

            scalarField& TVs       = TVibSpecies[speciei].primitiveFieldRef();
            const scalarField& eVs = eVibSpecies[speciei].primitiveField();
            const scalarField& rYs = rhoYi[speciei].primitiveField();

            forAll(TVs, celli)
            {
                if (rYs[celli] < vibHelpers::YTrace*rhoCells[celli]) continue;

                if (eVs[celli] <= VSMALL)
                {
                    TVs[celli] = TTRCells[celli];
                    ++nNonPos;
                    continue;
                }

                TVs[celli] =
                    vibHelpers::TvHO(eVs[celli], RSp[speciei], thetaSp[speciei]);

                if (debugLevel >= 2)
                {
                    const scalar eBack = mix.EsVib
                    (
                        speciei, pCells[celli], TTRCells[celli], TVs[celli],
                        thetaSp[speciei]
                    );
                    maxRT = max(maxRT, mag(eBack - eVs[celli])/eVs[celli]);
                }
            }
        }

        // (b) Cells: energy-equivalent mixture T_v over the molecules present
        forAll(TVibMixCells, celli)
        {
            label n = 0;
            forAll(Y_species, speciei)
            {
                if (!vibSpecieActive[speciei]) continue;

                const scalar Ys  = rhoYi[speciei][celli]/rhoCells[celli];
                const scalar evs = eVibSpecies[speciei][celli];
                if (Ys < vibHelpers::YTrace || evs <= VSMALL) continue;

                wY[n]  = Ys;
                wR[n]  = RSp[speciei];
                wTh[n] = thetaSp[speciei];
                wTv[n] = TVibSpecies[speciei][celli];
                wEv[n] = evs;
                ++n;
            }

            TVibMixCells[celli] = vibHelpers::TvMix
            (
                n, wY, wR, wTh, wTv, wEv, TTRCells[celli], maxG2
            );
        }

        // (c) Cells, trace molecules (G4): T_v = mixture T_v, e_v = e_v(T_v).
        //     eps_m (rhoYEvibSpecies) is not touched.
        forAll(Y_species, speciei)
        {
            if (!vibSpecieActive[speciei]) continue;

            scalarField& TVs       = TVibSpecies[speciei].primitiveFieldRef();
            scalarField& eVs       = eVibSpecies[speciei].primitiveFieldRef();
            const scalarField& rYs = rhoYi[speciei].primitiveField();

            forAll(TVs, celli)
            {
                if (rYs[celli] >= vibHelpers::YTrace*rhoCells[celli]) continue;

                TVs[celli] = TVibMixCells[celli];
                eVs[celli] = mix.EsVib
                (
                    speciei, pCells[celli], TTRCells[celli], TVs[celli],
                    thetaSp[speciei]
                );
                ++nTrace;
            }
        }

        // (d) Boundary faces of patches that do not fix T_v: same three steps.
        //     Fixed-value and jump (mixed) T_v patches are left to their BCs.
        forAll(TVib.boundaryField(), patchi)
        {
            fvPatchScalarField& pTVmix = TVib.boundaryFieldRef()[patchi];
            if (pTVmix.fixesValue()) continue;

            const fvPatchScalarField& pTTR = TTR.boundaryField()[patchi];

            forAll(pTVmix, facei)
            {
                label n = 0;
                forAll(Y_species, speciei)
                {
                    if (!vibSpecieActive[speciei]) continue;

                    fvPatchScalarField& pTVs =
                        TVibSpecies[speciei].boundaryFieldRef()[patchi];
                    const scalar Ys =
                        Y_species[speciei].boundaryField()[patchi][facei];
                    const scalar evs =
                        eVibSpecies[speciei].boundaryField()[patchi][facei];

                    if (Ys < vibHelpers::YTrace) continue;

                    // A molecule whose own BC fixes T_v on this
                    // patch (per-molecule jump, neSmoluchowskiJumpTv) keeps
                    // it; its energy for G2 is then e_v(T_v,face)
                    const bool fixedTv = pTVs.fixesValue();
                    scalar evUse = evs;

                    if (fixedTv)
                    {
                        evUse = vibHelpers::evHO
                        (
                            pTVs[facei], RSp[speciei], thetaSp[speciei]
                        );
                    }
                    else
                    {
                        if (evs <= VSMALL)
                        {
                            pTVs[facei] = pTTR[facei];
                            continue;
                        }

                        pTVs[facei] =
                            vibHelpers::TvHO(evs, RSp[speciei], thetaSp[speciei]);
                    }

                    if (evUse <= VSMALL) continue;

                    wY[n]  = Ys;
                    wR[n]  = RSp[speciei];
                    wTh[n] = thetaSp[speciei];
                    wTv[n] = pTVs[facei];
                    wEv[n] = evUse;
                    ++n;
                }

                pTVmix[facei] = vibHelpers::TvMix
                (
                    n, wY, wR, wTh, wTv, wEv, pTTR[facei], maxG2f
                );

                forAll(Y_species, speciei)
                {
                    if (!vibSpecieActive[speciei]) continue;
                    if
                    (
                        Y_species[speciei].boundaryField()[patchi][facei]
                    < vibHelpers::YTrace
                    && !TVibSpecies[speciei].boundaryField()[patchi].fixesValue()
                    )
                    {
                        TVibSpecies[speciei].boundaryFieldRef()[patchi][facei] =
                            pTVmix[facei];
                    }
                }
            }
        }

        forAll(Y_species, speciei)
        {
            if (!vibSpecieActive[speciei]) continue;
            TVibSpecies[speciei].correctBoundaryConditions();
            eVibSpecies[speciei].correctBoundaryConditions();
        }
        TVib.correctBoundaryConditions();
        if (debugLevel >= 2)
        {
            Info<< "  [vib] G1 round-trip max = "
                << returnReduce(maxRT, maxOp<scalar>())
                << "  G2 mismatch max = " << returnReduce(maxG2, maxOp<scalar>())
                << "  trace (G4) = " << returnReduce(nTrace, sumOp<label>())
                << "  e_v<=0 -> T_TR = " << returnReduce(nNonPos, sumOp<label>())
                << nl;
        }

        // ---- end P3-03 -------------------------------------------------------

        // Mutation++ path: atomic electronic energy diagnostic
        if (thermo2T.hasElectronicVibCoupling())
        {
            forAll(Y_species, speciei)
            {
                if (mix.isSpecieMolecular(speciei)) continue;

                scalarField& eVsAtom  = eVibSpecies[speciei].primitiveFieldRef();
                scalarField& rYeVAtom = rhoYEvibSpecies[speciei].primitiveFieldRef();
                const scalarField& TVibCells = TVib.primitiveField();

                forAll(TVibCells, celli)
                {
                    eVsAtom[celli] = mix.EsVib
                    (
                        speciei, p[celli], TTR[celli], TVibCells[celli], 0.0
                    );
                    rYeVAtom[celli] = rhoYi[speciei][celli]*eVsAtom[celli];
                }

                eVibSpecies[speciei].correctBoundaryConditions();
                rhoYEvibSpecies[speciei].correctBoundaryConditions();
            }

            // Recompute total eV including atomic electronic contribution
            eV = dimensionedScalar("zero", eV.dimensions(), 0.0);
            forAll(Y_species, speciei)
            {
                if (vibSpecieActive[speciei] || !mix.isSpecieMolecular(speciei))
                {
                    eV += rhoYEvibSpecies[speciei]/rho;
                }
            }
            eV.correctBoundaryConditions();
            rhoEv = rho*eV;
            rhoEv.correctBoundaryConditions();
        }
    }

    // ============================================================
    // Stage 3 — Viscous correction: vibrational conduction (G3)
    //
    // P3-08 (D3.1). Replaces the mixture form with the face-chord
    // diffusivity and the redistribution (K5, F2.13, F3.1).
    //
    // Per molecule m, Wilke-weighted with the same phi_m as kappa_TR:
    //     q_m = - X_m (kappa_v,m/phi_m) grad T_v,m
    //         = - X_m Gt_m grad e_v,m ,   Gt_m = (kappa_v,m/c_vv,m)/phi_m
    // (chain rule, exact for any kappa_v model). For the Eucken-type
    // kappa_v = mu c_vv/Sc of blottnerEucken2T, kappa_v/c_vv = mu/Sc: no
    // T_v dependence, no 0/0 in cold gas, no ratio of two gradients.
    //
    // Solved implicitly on y_m = eps_m/rho = Y_m e_v,m, using
    //     X_m Gt_m grad e_v,m = Gh_m (grad y_m - e_v,m grad Y_m),
    //     Gh_m = (W/W_m) Gt_m ,
    // the second part explicit. Nothing divides by Y_m or X_m.
    // The conducted energy is added to rhoE, as before.
    // ============================================================
    if (!inviscid)
    {
        const label nSp = Y_species.size();
        const scalarField& rhoC = rho.primitiveField();
        const scalarField& V    = mesh.V();

        scalarList WSp(nSp);
        scalarList thSp(nSp, Zero);
        forAll(Y_species, speciei)
        {
            WSp[speciei] = mix.W(speciei);
            if (vibSpecieActive[speciei]) thSp[speciei] = mix.thetaVib(speciei);
        }

        // Work lists for one point (cell or face)
        scalarList Xp(nSp), mup(nSp);

        // Gh_m at one point; Xp and mup must hold that point's values
        auto GhPoint = [&]
        (
            const label m,
            const scalar pP,
            const scalar TTRP,
            const scalar TvP,
            const scalar Wmix
        ) -> scalar
        {
            // Wilke denominator phi_m = sum_r X_r phi_mr (r = m gives X_m),
            // as WilkeMR::scalingFactor
            scalar phi = 0;
            for (label r = 0; r < nSp; ++r)
            {
                const scalar a =
                    1 + std::sqrt(mup[m]/mup[r])*std::pow(WSp[r]/WSp[m], 0.25);
                phi += Xp[r]*a*a/std::sqrt(8*(1 + WSp[m]/WSp[r]));
            }

            // kappa_v/c_vv at the local T_v; at T_v = theta where c_vv
            // underflows (frozen gas): the ratio is T_v-independent for Eucken
            scalar Tev = TvP;
            scalar cvv = mix.CvVib(m, pP, TTRP, Tev, thSp[m]);
            if (!(cvv > VSMALL))
            {
                Tev = thSp[m];
                cvv = mix.CvVib(m, pP, TTRP, Tev, thSp[m]);
            }
            const scalar kOverC =
                mix.kappaVib(m, pP, TTRP, Tev, thSp[m])/cvv;

            return (Wmix/WSp[m])*kOverC/max(phi, VSMALL);
        };

        // Fill Xp, mup at one point from mass fractions; returns W_mix
        auto fillPoint = [&]
        (
            const UList<scalar>& Ypt,
            const scalar TTRP
        ) -> scalar
        {
            scalar sumYW = 0;
            for (label r = 0; r < nSp; ++r) sumYW += Ypt[r]/WSp[r];
            const scalar Wmix = 1/max(sumYW, VSMALL);
            for (label r = 0; r < nSp; ++r)
            {
                Xp[r]  = Ypt[r]*Wmix/WSp[r];
                mup[r] = mix.mu(r, TTRP);
            }
            return Wmix;
        };

        scalarList Ypt(nSp);

        volScalarField dRhoE
        (
            IOobject("dRhoEvibCond", runTime.timeName(), mesh,
                    IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
            mesh,
            dimensionedScalar("zero", dimEnergy/dimVolume, 0.0)
        );

        autoPtr<volScalarField> eV_beforePtr;
        if (debugLevel >= 2)
        {
            eV_beforePtr.reset(new volScalarField("eV_before", eV));
        }

        scalar dEcond = 0;      // conducted vibrational energy, all molecules [J]

        // P4-07: vibrational wall heat flux, accumulated per molecule below
        forAll(qWallVib.boundaryField(), patchi)
        {
            qWallVib.boundaryFieldRef()[patchi] = 0.0;
        }

        forAll(Y_species, m)
        {
            if (!vibSpecieActive[m]) continue;

            // --- Gh_m: cells, processor patches, then physical patches
            volScalarField Gh
            (
                IOobject("Gh_" + Y_species[m].name(), runTime.timeName(), mesh,
                        IOobject::NO_READ, IOobject::NO_WRITE,
                        IOobject::NO_REGISTER),
                mesh,
                dimensionedScalar("zero", dimMass/dimLength/dimTime, 0.0)
            );
            {
                scalarField& GhC = Gh.primitiveFieldRef();
                const scalarField& pC   = p.primitiveField();
                const scalarField& TTRC = TTR.primitiveField();
                const scalarField& TvC  = TVibSpecies[m].primitiveField();

                forAll(GhC, celli)
                {
                    for (label r = 0; r < nSp; ++r)
                    {
                        Ypt[r] = Y_species[r][celli];
                    }
                    const scalar Wmix = fillPoint(Ypt, TTRC[celli]);
                    GhC[celli] =
                        GhPoint(m, pC[celli], TTRC[celli], TvC[celli], Wmix);
                }
            }
            Gh.correctBoundaryConditions();     // processor patches
            
            forAll(Gh.boundaryField(), patchi)
            {
                fvPatchScalarField& pGh = Gh.boundaryFieldRef()[patchi];
                if (pGh.coupled()) continue;

                const fvPatchScalarField& pp   = p.boundaryField()[patchi];
                const fvPatchScalarField& pTTR = TTR.boundaryField()[patchi];
                const fvPatchScalarField& pTv  =
                    TVibSpecies[m].boundaryField()[patchi];

                forAll(pGh, facei)
                {
                    for (label r = 0; r < nSp; ++r)
                    {
                        Ypt[r] = Y_species[r].boundaryField()[patchi][facei];
                    }
                    const scalar Wmix = fillPoint(Ypt, pTTR[facei]);
                    pGh[facei] =
                        GhPoint(m, pp[facei], pTTR[facei], pTv[facei], Wmix);
                }
            }

            // P4-06 (D4.4): Gamma_v,m = X_m G_m = Y_m Gh_m for the
            // per-molecule T_v jump BC (used at its next update)
            GammaVibSpecies[m] = Y_species[m]*Gh;

            // --- y_m = eps_m/rho with BCs: fixed where T_v is fixed (wall,
            //     inlet, T_v jump), zero gradient elsewhere; constraint and
            //     processor patches keep their own type
            const volScalarField y0("y0", rhoYEvibSpecies[m]/rho);

            wordList yTypes(mesh.boundary().size());
            forAll(yTypes, patchi)
            {
                yTypes[patchi] =
                    TVibSpecies[m].boundaryField()[patchi].fixesValue()
                ? fixedValueFvPatchScalarField::typeName
                : zeroGradientFvPatchScalarField::typeName;
            }

            volScalarField y
            (
                IOobject("eVib", runTime.timeName(), mesh,
                        IOobject::NO_READ, IOobject::NO_WRITE,
                        IOobject::NO_REGISTER),
                y0,
                yTypes
            );

            forAll(y.boundaryField(), patchi)
            {
                fvPatchScalarField& py = y.boundaryFieldRef()[patchi];
                if (!py.fixesValue()) continue;

                const fvPatchScalarField& pp   = p.boundaryField()[patchi];
                const fvPatchScalarField& pTTR = TTR.boundaryField()[patchi];
                const fvPatchScalarField& pTv  =
                    TVibSpecies[m].boundaryField()[patchi];
                const fvPatchScalarField& pY   =
                    Y_species[m].boundaryField()[patchi];

                scalarField yb(py.size());
                forAll(yb, facei)
                {
                    yb[facei] = pY[facei]*mix.EsVib
                    (
                        m, pp[facei], pTTR[facei], pTv[facei], thSp[m]
                    );
                }
                py == yb;
            }
            y.correctBoundaryConditions();

            // --- Implicit conduction (deferred-correction ddt form, as before)
            solve
            (
                fvm::ddt(rho, y) - fvc::ddt(rho, y)
            - fvm::laplacian(fvc::interpolate(Gh), y)
            ==
            - fvc::laplacian
                (
                    fvc::interpolate(Gh*eVibSpecies[m]),
                    Y_species[m]
                )
            );
            
            // P4-07 (D4.5): vibrational wall heat flux of this molecule from
            // the discrete flux of the solve above, positive into the wall
            {
                const surfaceScalarField GhF(fvc::interpolate(Gh));
                const surfaceScalarField GhEF
                (
                    fvc::interpolate(Gh*eVibSpecies[m])
                );
                const surfaceScalarField sgY(fvc::snGrad(y));
                const surfaceScalarField sgYm(fvc::snGrad(Y_species[m]));

                forAll(mesh.boundary(), patchi)
                {
                    if (!isA<wallFvPatch>(mesh.boundary()[patchi])) continue;

                    qWallVib.boundaryFieldRef()[patchi] +=
                        GhEF.boundaryField()[patchi]*sgYm.boundaryField()[patchi]
                    - GhF.boundaryField()[patchi]*sgY.boundaryField()[patchi];
                }
            }

            // --- Back to the conserved field; energy to rhoE
            const scalarField dEps(rhoC*(y.primitiveField() - y0.primitiveField()));

            rhoYEvibSpecies[m].primitiveFieldRef() += dEps;
            rhoYEvibSpecies[m].correctBoundaryConditions();

            eVibSpecies[m] =
                rhoYEvibSpecies[m]/max(rhoYi[m], vibHelpers::YTrace*rho);
            eVibSpecies[m].correctBoundaryConditions();

            dRhoE.primitiveFieldRef() += dEps;
            if (debugLevel >= 1)
            {
                dEcond += gSum(dEps*V);
            }

            if (debugLevel >= 2)
            {
                Info<< "  [G3] " << Y_species[m].name()
                    << "  Gh min/max = " << gMin(Gh.primitiveField())
                    << " / " << gMax(Gh.primitiveField()) << nl;
            }
        }

        // vibrational wall heat load [W]
        if (debugLevel >= 1)
        {
            scalar QV = 0;
            forAll(mesh.boundary(), patchi)
            {
                if (!isA<wallFvPatch>(mesh.boundary()[patchi])) continue;
                QV += sum
                (
                    qWallVib.boundaryField()[patchi]
                    *mesh.magSf().boundaryField()[patchi]
                );
            }
            Info<< "  [qWall] vib = " << returnReduce(QV, sumOp<scalar>())
                << " W" << nl;
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

        rhoE.primitiveFieldRef() += dRhoE.primitiveField();
        rhoE.correctBoundaryConditions();

        // ---------------------- diagnostics ----------------------
        if (debugLevel >= 2)
        {
            Info<< "  eV     min/max = " << gMin(eV.primitiveField())
                << " / " << gMax(eV.primitiveField()) << nl;

            forAll(Y_species, i)
            {
                if (!vibSpecieActive[i]) continue;
                Info<< "  " << Y_species[i].name()
                    << "  rhoYi min = " << gMin(rhoYi[i].primitiveField())
                    << "  eVib max = " << gMax(eVibSpecies[i].primitiveField())
                    << "  TVib max = " << gMax(TVibSpecies[i].primitiveField()) << nl;
            }
        }

        if (debugLevel >= 1)
        {
            const scalar Etot = gSum(rhoE.primitiveField()*V);
            const scalar Evib = gSum(rhoEv.primitiveField()*V);

            Info<< "  E_total    = " << Etot
                << "  E_vib      = " << Evib
                << "  E_vib/E_tot= " << Evib/Etot << nl
                << "  dE_vib(conduction) = " << dEcond
                << "   (should be small: boundary flux only)" << nl;
        }

        if (debugLevel >= 2)
        {
            Info<< "  max|delta_eV| = " << gMax(mag(eV.primitiveField()
                                            - eV_beforePtr().primitiveField())) << nl;
        }
    }
    // ---- end P3-08 -----------------------------------------------------------
}

// ************************************************************************* //
