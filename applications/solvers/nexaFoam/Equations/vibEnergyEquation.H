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

        PtrList<volScalarField> QCVempty;
        PtrList<volScalarField>* QCVptr = &QCVempty;
        if (useChemistry)
        {
            const PtrList<volScalarField::Internal>& RRchem = chemistryPtr->RR();
            QCVptr = &thermo2T.correctCVSource(RRchem);
        }
        PtrList<volScalarField>& QCV = *QCVptr;

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
            const volScalarField& QCVi = (speciei < QCV.size() ? QCV[speciei] : QVzero);

            solve
            (
                fvm::ddt(rhoYEvibSpecies[speciei])
              + fvc::div(phiRhoEvibYi[speciei] + phiJEvib_i)
             ==
                QVTi + QVVi + QCVi
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

                const scalar eBack = mix.EsVib
                (
                    speciei, pCells[celli], TTRCells[celli], TVs[celli],
                    thetaSp[speciei]
                );
                maxRT = max(maxRT, mag(eBack - eVs[celli])/eVs[celli]);
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
                    if (evs <= VSMALL)
                    {
                        pTVs[facei] = pTTR[facei];
                        continue;
                    }

                    pTVs[facei] =
                        vibHelpers::TvHO(evs, RSp[speciei], thetaSp[speciei]);

                    wY[n]  = Ys;
                    wR[n]  = RSp[speciei];
                    wTh[n] = thetaSp[speciei];
                    wTv[n] = pTVs[facei];
                    wEv[n] = evs;
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

        Info<< "  [vib] G1 round-trip max = "
            << returnReduce(maxRT, maxOp<scalar>())
            << "  G2 mismatch max = " << returnReduce(maxG2, maxOp<scalar>())
            << "  trace (G4) = " << returnReduce(nTrace, sumOp<label>())
            << "  e_v<=0 -> T_TR = " << returnReduce(nNonPos, sumOp<label>())
            << nl;
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
    // Stage 3 — Viscous correction (vibrational conduction)
    //
    // The physical term is  div(kappaVe * grad(TVib)).
    //
    // fvm::laplacian can only be implicit in the variable being solved, which
    // here is eV, so the term must be written  div(Gamma * grad(eV))  with
    //
    //     Gamma = kappaVe * dTVib/deV .
    //
    // The previous form used Gamma = kappaVe/CvVibMix with a CELL-CENTRED
    // CvVib. That is only equivalent when CvVib is smooth across a face.
    // CvVib = (R/W)(x/sinh x)^2, x = theta/(2 TVib), is exponential in the
    // near-wall range: across one wall-adjacent face in the expansion region
    // it goes from ~1 to ~152 J/kgK (145x) while the chord dTVib/deV across
    // that same face corresponds to ~69 J/kgK. Neither cell value is close,
    // and no interpolation of the two recovers it.
    //
    // Instead build Gamma directly at the face from the two gradients:
    //
    //     Gamma_f = kappaVe_f * snGrad(TVib) / snGrad(eV)
    //
    // so that  Gamma_f*snGrad(eV) == kappaVe_f*snGrad(TVib)  identically, for
    // any variation of CvVib. This is exact, fully implicit, and positive by
    // construction (eV is monotone in TVib), so it carries no diffusive
    // time-step penalty. CvVibMix is retained only as the fallback where both
    // gradients vanish, and for diagnostics.
    //
    // The equivalent substitution is harmless for eTR because CvTR ~ 2.5R is
    // very nearly constant; it is not harmless here.
    // ============================================================
    if (!inviscid)
    {
        // --- Mixture vibrational specific heat (fallback + diagnostic) ---
        volScalarField CvVibMix
        (
            IOobject("CvVibMix", runTime.timeName(), mesh,
                    IOobject::NO_READ, IOobject::AUTO_WRITE, IOobject::NO_REGISTER),
            mesh,
            dimensionedScalar("small", dimEnergy/dimMass/dimTemperature, SMALL),
            calculatedFvPatchScalarField::typeName
        );
        {
            scalarField& cvv = CvVibMix.primitiveFieldRef();
            forAll(cvv, celli)
            {
                scalar sumc = 0.0;
                forAll(Y_species, speciei)
                {
                    if (!vibSpecieActive[speciei])       continue;
                    if (!mix.isSpecieMolecular(speciei)) continue;   // atomic: no vib mode

                    const scalar Tv = max(TVibSpecies[speciei][celli], scalar(1));
                    const scalar x  = 0.5*mix.thetaVib(speciei)/Tv;
                    const scalar xs = (x < scalar(20)) ? x/Foam::sinh(x) : 0.0;  // ->0 as Tv->0, no sinh overflow

                    // CvVib_s = (R/M_s)*(x/sinh x)^2  — identical to the thermo's SHO CvVib
                    sumc += Y_species[speciei][celli]*mix.R(speciei)*xs*xs;
                }
                cvv[celli] = max(sumc, SMALL);
            }
            volScalarField::Boundary& cvvBf = CvVibMix.boundaryFieldRef();

            forAll(cvvBf, patchi)
            {
                fvPatchScalarField& cvvp = cvvBf[patchi];

                forAll(cvvp, facei)
                {
                    scalar sumc = 0.0;

                    forAll(Y_species, speciei)
                    {
                        if (!vibSpecieActive[speciei])       continue;
                        if (!mix.isSpecieMolecular(speciei)) continue;

                        const scalar Tv = max
                        (
                            TVibSpecies[speciei].boundaryField()[patchi][facei],
                            scalar(1)
                        );
                        const scalar x  = 0.5*mix.thetaVib(speciei)/Tv;
                        const scalar xs = (x < scalar(20)) ? x/Foam::sinh(x) : 0.0;

                        const scalar Ys = max
                        (
                            Y_species[speciei].boundaryField()[patchi][facei],
                            scalar(0)
                        );

                        sumc += Ys*mix.R(speciei)*xs*xs;
                    }

                    cvvp[facei] = max(sumc, SMALL);
                }
            }
        }

        volScalarField kappaVibMix("kappaVibMix", thermo2T.kappaVib());

        // Cell-based diffusivity: fallback only, where the face gradients are
        // too small for the chord to be meaningful.
        surfaceScalarField alphaVibF
        (
            "alphaVibF",
            fvc::interpolate(kappaVibMix/CvVibMix)
        );

        // --- Face-chord diffusivity: Gamma_f = kappaVe_f * dTVib_f / deV_f ---
        {
            const surfaceScalarField kappaVeF(fvc::interpolate(kappaVibMix));
            const surfaceScalarField sgT(fvc::snGrad(TVib));
            const surfaceScalarField sgE(fvc::snGrad(eV));

            // Below these the chord is numerically meaningless; keep fallback.
            const scalar sgTMin = 1e-6;    // [K/m]
            const scalar sgEMin = 1e-6;    // [J/kg/m]

            scalarField& af = alphaVibF.primitiveFieldRef();
            const scalarField& kf = kappaVeF.primitiveField();
            const scalarField& gt = sgT.primitiveField();
            const scalarField& ge = sgE.primitiveField();

            forAll(af, facei)
            {
                if (mag(gt[facei]) > sgTMin && mag(ge[facei]) > sgEMin)
                {
                    // mag() on both: eV is monotone in TVib, so the two
                    // gradients always share a sign and Gamma_f > 0.
                    af[facei] = kf[facei]*mag(gt[facei])/mag(ge[facei]);
                }
            }

            surfaceScalarField::Boundary& afBf = alphaVibF.boundaryFieldRef();
            const surfaceScalarField::Boundary& kfBf = kappaVeF.boundaryField();
            const surfaceScalarField::Boundary& gtBf = sgT.boundaryField();
            const surfaceScalarField::Boundary& geBf = sgE.boundaryField();

            forAll(afBf, patchi)
            {
                fvsPatchScalarField& afp = afBf[patchi];

                forAll(afp, facei)
                {
                    const scalar gtp = gtBf[patchi][facei];
                    const scalar gep = geBf[patchi][facei];

                    if (mag(gtp) > sgTMin && mag(gep) > sgEMin)
                    {
                        afp[facei] = kfBf[patchi][facei]*mag(gtp)/mag(gep);
                    }
                }
            }
        }

        // Capture eV before the conduction correction (for per-species
        // redistribution and for the energy diagnostics below)
        volScalarField eV_before("eV_before", eV);

        // Implicit vibrational conduction, deferred-correction ddt form
        solve
        (
            fvm::ddt(rho, eV) - fvc::ddt(rho, eV)
          - fvm::laplacian(alphaVibF, eV)
        );

        eV.correctBoundaryConditions();
        rhoEv = rho*eV;
        rhoEv.correctBoundaryConditions();

        rhoE += rho*(eV - eV_before);

        if (mixtureCheck)
        {
            volScalarField delta_eV("delta_eV", eV - eV_before);

            volScalarField Yactive
            (
                IOobject("Yactive", runTime.timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE,
                         IOobject::NO_REGISTER),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            );

            forAll(Y_species, speciei)
            {
                if (!vibSpecieActive[speciei]) continue;
                Yactive += Y_species[speciei];
            }

            const dimensionedScalar YactiveMin("small", dimless, SMALL);

            forAll(Y_species, speciei)
            {
                if (!vibSpecieActive[speciei]) continue;

                rhoYEvibSpecies[speciei] +=
                    rhoYi[speciei]*delta_eV/max(Yactive, YactiveMin);

                rhoYEvibSpecies[speciei].correctBoundaryConditions();

                eVibSpecies[speciei] =
                    rhoYEvibSpecies[speciei]
                    /max(rhoYi[speciei], vibHelpers::YTrace*rho);   // P3-02: no rhoYi floor
                eVibSpecies[speciei].correctBoundaryConditions();
            }

            // Round-trip check: this should now be ~1.0 to machine precision.
            volScalarField eVrebuilt
            (
                IOobject("eVrebuilt", runTime.timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE,
                         IOobject::NO_REGISTER),
                mesh,
                dimensionedScalar("zero", eV.dimensions(), 0.0)
            );

            forAll(Y_species, speciei)
            {
                if (!vibSpecieActive[speciei]) continue;
                eVrebuilt += rhoYEvibSpecies[speciei]/rho;
            }

            const scalar dNew =
                gSum((rho*(eVrebuilt - eV_before))().primitiveField()*mesh.V());
            const scalar dRef =
                gSum((rho*delta_eV)().primitiveField()*mesh.V());

            Info<< "  [vibRedist] recovered/applied = "
                << (mag(dRef) > VSMALL ? dNew/dRef : 1.0) << nl;
        }

        // ---------------------- diagnostics ----------------------
        Info<< "  eV     min/max = " << gMin(eV.primitiveField())
            << " / " << gMax(eV.primitiveField()) << nl;
        Info<< "  CvVibMix min   = " << gMin(CvVibMix.primitiveField()) << nl;
        Info<< "  alphaVibF max  = " << gMax(alphaVibF.primitiveField())
            << "   (cell-based would give "
            << gMax((kappaVibMix/CvVibMix)().primitiveField()) << ")" << nl;

        forAll(Y_species, i)
        {
            if (!vibSpecieActive[i]) continue;
            Info<< "  " << Y_species[i].name()
                << "  rhoYi min = " << gMin(rhoYi[i].primitiveField())
                << "  eVib max = " << gMax(eVibSpecies[i].primitiveField())
                << "  TVib max = " << gMax(TVibSpecies[i].primitiveField()) << nl;
        }

        const scalarField& V = mesh.V();
        const scalar Etot = gSum(rhoE.primitiveField()*V);
        const scalar Evib = gSum(rhoEv.primitiveField()*V);
        const scalar dEvibCond = gSum(rho.primitiveField()
                               *(eV.primitiveField() - eV_before.primitiveField())*V);

        Info<< "  E_total    = " << Etot
            << "  E_vib      = " << Evib
            << "  E_vib/E_tot= " << Evib/Etot << nl
            << "  dE_vib(conduction) = " << dEvibCond
            << "   (should be small: boundary flux only)" << nl;

        Info<< "  max|delta_eV| = " << gMax(mag(eV.primitiveField()
                                        - eV_before.primitiveField())) << nl;
    }
}

// ************************************************************************* //
