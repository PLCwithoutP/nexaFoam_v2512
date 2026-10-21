if (mixtureCheck)
{
    PtrList<volScalarField>& Y = thermo2T.composition().Y();
    const PtrList<surfaceScalarField>& phiRhoYi = kurganovFlux.phiRhoYi();

    const label nSpec = Y.size();
    const label last  = nSpec - 1;

    // --- 1) Solve partial densities rhoYi  (CONVECTION ONLY, explicit KNP)
    //     Mass diffusion removed from here; handled implicitly in step 4.
    for (label i = 0; i < nSpec; ++i)
    {
        solve
        (
            fvm::ddt(rhoYi[i])
          + fvc::div(phiRhoYi[i])
        );
    }

    // --- 2) Build mixture rho from the convected rhoYi (internal cells).
    //     Mass diffusion redistributes species but does NOT change mixture
    //     rho, so rho is frozen for the corrector below.
    rho.primitiveFieldRef() = scalar(0);
    for (label i = 0; i < nSpec; ++i)
    {
        rho.primitiveFieldRef() += rhoYi[i].primitiveField();
    }
    rho.primitiveFieldRef() =
        max(rho.primitiveField(), scalarField(rho.size(), SMALL));

    // --- 3) Recover Yi = rhoYi / rho  (no clipping yet; let diffusion act
    //     on the raw convected field, bound afterwards in step 5).
    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] = rhoYi[i]/rho;
        Y[i].correctBoundaryConditions();
    }

    // --- 4) IMPLICIT mass-diffusion corrector on Yi.
    //     Same deferred-correction form as the eTR / eV conduction solves.
    //     Shared mixture Deff => sum of the species Laplacians is the
    //     Laplacian of (sum Yi)=1 = 0, so all nSpec may be solved and the
    //     exact closure is restored by the renormalisation in step 5.
    for (label i = 0; i < nSpec; ++i)
    {
        solve
        (
            fvm::ddt(rho, Y[i]) - fvc::ddt(rho, Y[i])
          - fvm::laplacian(rho*Deff, Y[i])
        );
    }

    // --- 5) Conserved partial densities from the diffused Yi (P3-04, D3.2).
    //     The diffused Yi are renormalised (this only redistributes the linear
    //     solver residual; sum rhoYi = rho exactly) but NOT clipped: an
    //     undershoot stays in rhoYi so that every species mass is conserved.
    //     It is counted, and clipped only in the property field Y (step 6).
    volScalarField sumY
    (
        IOobject("sumY", runTime.timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar(dimless, 0.0)
    );
    for (label i = 0; i < nSpec; ++i) sumY += Y[i];
    sumY.max(SMALL);

    label  nNegY   = 0;
    scalar negMass = 0;
    for (label i = 0; i < nSpec; ++i)
    {
        rhoYi[i] = rho*Y[i]/sumY;
        rhoYi[i].correctBoundaryConditions();

        const scalarField& rYi = rhoYi[i].primitiveField();
        forAll(rYi, celli)
        {
            if (rYi[celli] < -1e-12*rho[celli])
            {
                ++nNegY;
                negMass -= rYi[celli]*mesh.V()[celli];
            }
        }
    }
    reduce(nNegY, sumOp<label>());
    reduce(negMass, sumOp<scalar>());
    if (nNegY > 0)
    {
        Info<< "    [Y] rhoY_s < -1e-12 rho kept in the conserved field: " << nNegY
            << " species-cells, " << negMass << " kg in total" << nl;
    }

    // --- 6) Property mass fractions (read by thermo, transport, chemistry and
    //     the fluxes): clip to [0, 1] and renormalise. rhoYi is not touched.
    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] = rhoYi[i]/rho;
        Y[i].max(0.0);
        Y[i].min(1.0);
    }
    sumY = dimensionedScalar(dimless, 0.0);
    for (label i = 0; i < nSpec; ++i) sumY += Y[i];
    sumY.max(SMALL);
    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] /= sumY;
        Y[i].correctBoundaryConditions();
    }
    // ---- end P3-04 ----------------------------------------------------------

    // --- 7) Rebuild diffusive mass flux phiJYi from the CORRECTED Yi for the
    //     explicit species-enthalpy / vib-energy transport in the energy
    //     equations (those terms remain explicit by design). N-1 from Fick,
    //     last closed as -sum so that sum(phiJYi) = 0 exactly.
    forAll(phiJYi, i)
    {
        phiJYi[i] =
            dimensionedScalar("zero", phiJYi[i].dimensions(), 0.0);
    }

    surfaceScalarField sumPhiJ
    (
        IOobject("sumPhiJ", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar(dimMass/dimTime, 0.0)
    );
    sumPhiJ.setOriented(true);

    for (label i = 0; i < last; ++i)
    {
        phiJYi[i] = Foam::FickModel::phiJ(rho, Deff, Y[i]);
        sumPhiJ  += phiJYi[i];
    }
    phiJYi[last] = -sumPhiJ;
}
