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
                const scalar eVsClamped = max(eVs[celli], scalar(0));

                if (eVsClamped < SMALL)
                {
                    TVs[celli] = scalar(0);
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
                const fvPatchScalarField& pp   = pBf[patchi];
                const fvPatchScalarField& pTTR = TTRBf[patchi];
                const fvPatchScalarField& peV  = eVBf[patchi];
                fvPatchScalarField& pTV        = TVBf[patchi];

                forAll(pp, facei)
                {
                    const scalar eVsClamped = max(peV[facei], scalar(0));

                    if (eVsClamped < SMALL)
                    {
                        pTV[facei] = scalar(0);
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

        forAll(TVib.boundaryField(), patchi)
        {
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
    // ============================================================
    if (!inviscid)
    {
        volScalarField kappaVibField("kappaVibField", thermo2T.kappaVib());

        surfaceScalarField qVibCorr
        (
            - fvc::interpolate(kappaVibField)
            * mesh.magSf()
            * fvc::snGrad(TVib)
        );
        qVibCorr.setOriented(true);

        // Capture eV before correction — needed for delta redistribution
        volScalarField eV_before("eV_before", eV);

        rhoEv -= runTime.deltaT() * fvc::div(qVibCorr);
        rhoEv.correctBoundaryConditions();
        eV = rhoEv/rho;
        eV.correctBoundaryConditions();

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
    }
}
