
// ---------------------------------------------------------------------------
// transportUpdate.H — nexaFoam Review Project, Phase 4 (P4-03; D4.1, D4.3)
//
// One Wilke-mixed kappa_TR for the TR conduction, kappaPP, the TR jump BC
// (registered field kappaTRW) and the Fick coefficient
//     D_eff = Le kappa_TR,W/(rho c_p,TR),   Le from thermophysicalProperties
// (C1, C3, C5, C6, K6). Wilke denominator phi_s = sum_r X_r phi_sr as
// WilkeMR::scalingFactor. Order (C10): cells -> processor patches ->
// physical faces from the face composition.
// Pure mixtures keep the library values.
// ---------------------------------------------------------------------------
if (mixtureCheck)
{
    const label nSp = Y_species.size();

    scalarList Ypt(nSp), Xpt(nSp), mupt(nSp), kpt(nSp), WSp(nSp), RSp(nSp);
    forAll(Y_species, s)
    {
        WSp[s] = mix.W(s);
        RSp[s] = mix.R(s);
    }

    // Wilke kappa_TR at one point from Ypt; returns kappa, sets cp and cv
    auto wilkePoint = [&]
    (
        const scalar pP,
        const scalar TP,
        scalar& cpMix,
        scalar& cvMix
    ) -> scalar
    {
        scalar sumYW = 0;
        for (label r = 0; r < nSp; ++r) sumYW += Ypt[r]/WSp[r];
        const scalar Wmix = 1/max(sumYW, VSMALL);

        cpMix = 0;
        cvMix = 0;
        for (label r = 0; r < nSp; ++r)
        {
            Xpt[r]  = Ypt[r]*Wmix/WSp[r];
            mupt[r] = mix.mu(r, TP);
            kpt[r]  = mix.kappaTR(r, pP, TP);
            const scalar cps = mix.CpTR(r, pP, TP);
            cpMix += Ypt[r]*cps;
            cvMix += Ypt[r]*(cps - RSp[r]);
        }

        scalar k = 0;
        for (label s = 0; s < nSp; ++s)
        {
            scalar phi = 0;
            for (label r = 0; r < nSp; ++r)
            {
                const scalar a =
                    1 + std::sqrt(mupt[s]/mupt[r])*std::pow(WSp[r]/WSp[s], 0.25);
                phi += Xpt[r]*a*a/std::sqrt(8*(1 + WSp[s]/WSp[r]));
            }
            k += Xpt[s]*kpt[s]/max(phi, VSMALL);
        }
        return k;
    };

    // Cells
    {
        const scalarField& pC   = p.primitiveField();
        const scalarField& TC   = TTR.primitiveField();
        const scalarField& rhoC = rho.primitiveField();
        scalarField& kC  = kappaTRW.primitiveFieldRef();
        scalarField& cvC = cvTRmix.primitiveFieldRef();
        scalarField& DC  = Deff.primitiveFieldRef();

        forAll(kC, celli)
        {
            for (label r = 0; r < nSp; ++r) Ypt[r] = Y_species[r][celli];
            scalar cp = 0, cv = 0;
            kC[celli]  = wilkePoint(pC[celli], TC[celli], cp, cv);
            cvC[celli] = cv;
            DC[celli]  = LeNum*kC[celli]/max(rhoC[celli]*cp, VSMALL);
        }
    }
    kappaTRW.correctBoundaryConditions();      // processor patches
    cvTRmix.correctBoundaryConditions();
    Deff.correctBoundaryConditions();

    // Physical faces from the face composition
    forAll(kappaTRW.boundaryField(), patchi)
    {
        fvPatchScalarField& pk = kappaTRW.boundaryFieldRef()[patchi];
        if (pk.coupled()) continue;

        fvPatchScalarField& pcv = cvTRmix.boundaryFieldRef()[patchi];
        fvPatchScalarField& pD  = Deff.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pp   = p.boundaryField()[patchi];
        const fvPatchScalarField& pT   = TTR.boundaryField()[patchi];
        const fvPatchScalarField& prho = rho.boundaryField()[patchi];

        forAll(pk, facei)
        {
            for (label r = 0; r < nSp; ++r)
            {
                Ypt[r] = Y_species[r].boundaryField()[patchi][facei];
            }
            scalar cp = 0, cv = 0;
            pk[facei]  = wilkePoint(pp[facei], pT[facei], cp, cv);
            pcv[facei] = cv;
            pD[facei]  = LeNum*pk[facei]/max(prho[facei]*cp, VSMALL);
        }
    }
}
else
{
    kappaTRW = thermo2T.kappaTR();
    Deff = thermo2T.fickDiffusionCoeff();
    Deff.correctBoundaryConditions();
}
// ---- end P4-03 -------------------------------------------------------------