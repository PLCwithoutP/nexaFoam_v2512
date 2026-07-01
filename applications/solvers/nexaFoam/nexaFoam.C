/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    nexaFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver based on
    central-upwind schemes of Kurganov and Tadmor with
    support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "fixed2TRhoFvPatchScalarField.H"
#include "maxwellSlipU.H"
#include "smoluchowskiJumpT.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "FluxCalculators/kurganov.H"
#include "DiffusionModels/FickModel.H"
#include "tcLibraryInterface.H"
#include "chemistryInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //#include "logDebugCreate.H"
    argList::addNote
    (
        "Density-based compressible flow solver based on"
        " central-upwind schemes of Kurganov and Tadmor with"
        " support for mesh-motion and topology changes."
    );

    #define NO_CONTROL
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    const dimensionedScalar v_zero(dimVolume/dimTime, Zero);

    // Flux Calculator
    autoPtr<Foam::kurganovFluxCalculator> kurganovFluxPtr(nullptr);

    if (!use2T)
    {
        kurganovFluxPtr.reset
        (
            new Foam::kurganovFluxCalculator
            (
                mesh,
                thermo2T,
                rho,
                rhoU,
                U,
                psi,
                eTR,
                TTR,
                pos,
                neg,
                phi,
                fluxScheme,
                v_zero,
                mixtureCheck,
                use2T
            )
        );
    }
    else if (use2T && !mixtureCheck)
    {
        kurganovFluxPtr.reset
        (
            new Foam::kurganovFluxCalculator
            (
                mesh,
                thermo2T,
                rho,
                rhoU,
                U,
                psi,
                eTR,
                TTR,
                eV,
                pos,
                neg,
                phi,
                fluxScheme,
                v_zero,
                mixtureCheck,
                use2T
            )
        );
    }
    else
    {
        kurganovFluxPtr.reset
        (
            new Foam::kurganovFluxCalculator
            (
                mesh,
                thermo2T,
                rho,
                rhoU,
                U,
                psi,
                eTR,
                TTR,
                eVibSpecies,
                pos,
                neg,
                phi,
                fluxScheme,
                v_zero,
                mixtureCheck,
                use2T
            )
        );
    }

    Foam::kurganovFluxCalculator& kurganovFlux = kurganovFluxPtr();

    surfaceVectorField& phiUp = kurganovFlux.phiUp();
    surfaceScalarField& phiEp = kurganovFlux.phiEp();
    surfaceScalarField& max_a = kurganovFlux.maxA();
    surfaceScalarField& phiEv = kurganovFlux.phiEv();
    const PtrList<surfaceScalarField>& phiRhoEvibYi = kurganovFlux.phiRhoEvibYi();

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (!LTS)
        {
            #include "setDeltaT.H"

            ++runTime;

            // Do any mesh changes
            mesh.update();
        }

        // Effective diffusion coefficient must be updated in every timestep
        Deff = thermo2T.fickDiffusionCoeff();
        Deff.correctBoundaryConditions();

        kurganovFlux.calcDirectedInterpolation();
        kurganovFlux.calcSoundSpeed();
        kurganovFlux.calcPropagationSpeedsAndWeights();
        kurganovFlux.calcAphi();

        #include "centralCourantNo.H"
        

        if (LTS)
        {
            #include "setRDeltaT.H"

            ++runTime;
        }

        
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (mixtureCheck)
        {
            kurganovFlux.calcSpeciesFluxes();
        }
        else
        {
            kurganovFlux.calcConservativeFluxes();
        }

        if (use2T)
        {
            if (mixtureCheck)
            {
                kurganovFlux.calcSpeciesVibFluxes();
            }
            else
            {
                kurganovFlux.calcPureMixtureVibFlux();
            }
        }

        muLam = thermo2T.mu();
        kappaPP = thermo2T.kappaTR();
        volTensorField tauMC("tauMC", muLam*dev2(Foam::T(fvc::grad(U))));

        if (use2T)
        {
            kappaVe = thermo2T.kappaVib();
        }

        #include "Equations/continuityEquation.H"

        #include "Equations/yEquation.H"

        #include "Equations/momentumEquation.H"

        #include "Equations/totalEnergyEquation.H"

        #include "Equations/vibEnergyEquation.H"

        // Sanity check
        if (mixtureCheck)
        {
            volScalarField drho(rhoCont - rho);
            const scalar maxErr = gMax(mag(drho)().primitiveField());
            Info<< "Sanity check: max(|rhoCont - rhoFromSpecies|) = "
                << maxErr << nl << endl;
        }

        #include "Equations/applyChemistry.H"

        // --- DEBUG: cylinder wall face flux (remove after diagnosis) ---
{
    const label patchi = mesh.boundaryMesh().findPatchID("cylinder");

    if (patchi != -1)
    {
        const fvsPatchScalarField& phiw = phi.boundaryField()[patchi];
        const fvPatchVectorField&  Uw   = U.boundaryField()[patchi];
        const fvPatchScalarField&  rhow = rho.boundaryField()[patchi];
        const fvPatchScalarField&  pw   = thermo2T.p().boundaryField()[patchi];
        const vectorField          Sf   = mesh.Sf().boundaryField()[patchi];
        const scalarField          magSf(mesh.magSf().boundaryField()[patchi]);
        const vectorField          nf(mesh.boundary()[patchi].nf());

        // Un-normalised (sum) and worst-face diagnostics
        const scalar sumPhi   = gSum(phiw);
        const scalar sumMagPhi= gSum(mag(phiw));
        const scalar maxPhi   = gMax(mag(phiw));

        Info<< "WALL cylinder:"
            << "  sum(phi)="    << sumPhi
            << "  sum|phi|="    << sumMagPhi
            << "  max|phi|="    << maxPhi << nl;

        // Per-face look at the first few faces (nose region)
        forAll(phiw, facei)
        {
            if (facei < 1500)
            {
                Info<< "   face" << facei
                    << "  phi="   << phiw[facei]
                    << "  U.n="   << (Uw[facei] & nf[facei])
                    << "  |U|="   << mag(Uw[facei])
                    << "  rho="   << rhow[facei]
                    << "  p="     << pw[facei] << nl;
            }
        }
        Info<< endl;
    }
}
// --------------------------------------------------------------

        runTime.write();

        runTime.printExecutionTime(Info);

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
