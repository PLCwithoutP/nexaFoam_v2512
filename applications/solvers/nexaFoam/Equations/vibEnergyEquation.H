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

            rhoYi[speciei].max(SMALL);
            rhoYEvibSpecies[speciei].primitiveFieldRef() =
                max(rhoYEvibSpecies[speciei].primitiveField(),
                    scalarField(rhoYEvibSpecies[speciei].size(), scalar(0)));
            rhoYEvibSpecies[speciei].correctBoundaryConditions();

            eVibSpecies[speciei] = rhoYEvibSpecies[speciei]/rhoYi[speciei];
            eVibSpecies[speciei].correctBoundaryConditions();
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
        const scalarField& pCells   = p.primitiveField();
        const scalarField& TTRCells = TTR.primitiveField();
        scalarField& TVibMixCells   = TVib.primitiveFieldRef();

        // Per-species TVib inversion
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

            forAll(TVs, celli)
            {
                const scalar eVsClamped = max(eVs[celli], scalar(20));

                if (eVsClamped < SMALL)
                {
                    TVs[celli] = TTRCells[celli];
                    continue;
                }

                TVs[celli] = mix.TEVib
                (
                    speciei,
                    eVsClamped,
                    pCells[celli],
                    TTRCells[celli],
                    TVs[celli]
                );
            }

            const volScalarField::Boundary& pBf   = p.boundaryField();
            const volScalarField::Boundary& TTRBf = TTR.boundaryField();
            const volScalarField::Boundary& eVBf  =
                eVibSpecies[speciei].boundaryField();
            volScalarField::Boundary& TVBf =
                TVibSpecies[speciei].boundaryFieldRef();

            forAll(pBf, patchi)
            {
                if (TVBf[patchi].fixesValue()) continue;
                const fvPatchScalarField& pp   = pBf[patchi];
                const fvPatchScalarField& pTTR = TTRBf[patchi];
                const fvPatchScalarField& peV  = eVBf[patchi];
                fvPatchScalarField& pTV        = TVBf[patchi];

                forAll(pp, facei)
                {
                    const scalar eVsClamped = max(peV[facei], scalar(20));

                    if (eVsClamped < SMALL)
                    {
                        pTV[facei] = pTTR[facei];
                        continue;
                    }

                    pTV[facei] = mix.TEVib
                    (
                        speciei,
                        eVsClamped,
                        pp[facei],
                        pTTR[facei],
                        pTV[facei]
                    );
                }
            }

            TVibSpecies[speciei].correctBoundaryConditions();
        }

        // Mixture TVib: mass-fraction-weighted average
        forAll(TVibMixCells, celli)
        {
            scalar sumY  = 0.0;
            scalar sumTV = 0.0;

            forAll(Y_species, speciei)
            {
                if (!vibSpecieActive[speciei]) continue;
                const scalar Yi = max(Y_species[speciei][celli], scalar(0));
                sumY  += Yi;
                sumTV += Yi*TVibSpecies[speciei][celli];
            }

            TVibMixCells[celli] = (sumY > SMALL ? sumTV/sumY : 0.0);
        }

        // NOTE: this writes the TVib patch values directly. On a fixedValue
        // patch the subsequent correctBoundaryConditions() does NOT restore
        // the dictionary value, so the wall TVib is whatever the per-species
        // inversion produced. Intentional for jump/slip walls; for a plain
        // fixedValue wall it silently overrides 0/TVib. See known issues.
        forAll(TVib.boundaryField(), patchi)
        {
            if (TVib.boundaryField()[patchi].fixesValue()) continue;

            fvPatchScalarField& pTVmix = TVib.boundaryFieldRef()[patchi];

            forAll(pTVmix, facei)
            {
                scalar sumY  = 0.0;
                scalar sumTV = 0.0;

                forAll(Y_species, speciei)
                {
                    if (!vibSpecieActive[speciei]) continue;

                    const scalar Yi = max(
                        Y_species[speciei].boundaryField()[patchi][facei],
                        scalar(0));
                    sumY  += Yi;
                    sumTV += Yi*TVibSpecies[speciei].boundaryField()[patchi][facei];
                }

                pTVmix[facei] = (sumY > SMALL ? sumTV/sumY : 0.0);
            }
        }

        TVib.correctBoundaryConditions();

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

        // KNOWN ISSUE: rhoE is not updated for the energy this term moves
        // between cells, so eTR = rhoE/rho - KE - eV converts it in place
        // instead of transporting it. See known issues; enabling this alone
        // was previously unstable and is still open.
        //rhoE += rho*(eV - eV_before);
        //rhoE.correctBoundaryConditions();

        if (mixtureCheck)
        {
            // Delta eV from vibrational conduction correction
            volScalarField delta_eV("delta_eV", eV - eV_before);

            forAll(Y_species, speciei)
            {
                if (!vibSpecieActive[speciei]) continue;

                rhoYEvibSpecies[speciei] += rhoYi[speciei]*delta_eV;
                rhoYEvibSpecies[speciei].correctBoundaryConditions();

                eVibSpecies[speciei] =
                    rhoYEvibSpecies[speciei]/rhoYi[speciei];
                eVibSpecies[speciei].correctBoundaryConditions();
            }
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
