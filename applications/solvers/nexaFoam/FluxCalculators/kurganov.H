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

Description
    Kurganov-type (central-upwind) face flux calculator.
    - Kurganov-Tadmor ("Tadmor") and Kurganov-Noelle-Petrova (default) variants
    - Mesh motion supported (relative flux)
    - Designed as a black-box: inputs/outputs are explicit

Notes
    This header is intended to be included by the solver translation unit.
    It defines a header-only class with inline functions.

\*---------------------------------------------------------------------------*/

#ifndef nexaFoam_KURGANOV_FLUX_CALCULATOR_H
#define nexaFoam_KURGANOV_FLUX_CALCULATOR_H

#include "fvCFD.H"
#include "directionInterpolate.H"
#include "tcLibraryInterface.H"
#include "PtrList.H"

namespace Foam
{

class kurganovFluxCalculator
{
    // --- Required references (inputs)
    const fvMesh& mesh_;
    tcLibraryInterface& thermo2T_;

    const volScalarField& rho_;
    const volVectorField& rhoU_;
    const volVectorField& U_;

    const volScalarField& psi_;
    const volScalarField& eTR_;
    const volScalarField& TTR_;

    const surfaceScalarField& pos_;
    const surfaceScalarField& neg_;


    surfaceScalarField& phi_;          // mixture mass flux (updated)

    const word& fluxScheme_;           // "Tadmor" or default KNP
    const dimensionedScalar& v_zero_;  // zero with dimVolume/dimTime

    // --- Intermediate face fields (workspace)
    surfaceScalarField rho_pos_;
    surfaceScalarField rho_neg_;

    surfaceVectorField rhoU_pos_;
    surfaceVectorField rhoU_neg_;

    surfaceScalarField rPsi_pos_;
    surfaceScalarField rPsi_neg_;

    surfaceScalarField e_pos_;
    surfaceScalarField e_neg_;

    surfaceVectorField U_pos_;
    surfaceVectorField U_neg_;

    surfaceScalarField p_pos_;
    surfaceScalarField p_neg_;

    surfaceScalarField phi_fp_;  // unoriented
    surfaceScalarField phi_fn_;  // unoriented

    surfaceScalarField cSf_pos_;
    surfaceScalarField cSf_neg_;

    surfaceScalarField a_fp_;
    surfaceScalarField a_fn_;

    surfaceScalarField alpha_;
    surfaceScalarField minus_alpha_;

    surfaceScalarField nw_f_;   // diffusive flux factor (negative for Tadmor)
    surfaceScalarField aphi_fp_;
    surfaceScalarField aphi_fn_;

    bool mixtureCheck_; // controls whether single specie or multi specie flux will be employed
    // --- Multispecies storage (only used when mixtureCheck_ == true)
    bool use2T_;

    label nSpecie_;
    PtrList<surfaceScalarField> phiRhoYi_;   // face convective fluxes of rho*Yi

    // Optional vibrational inputs
    const volScalarField* eVPtr_;                    // pure-mixture 2T
    const PtrList<volScalarField>* eVibSpeciesPtr_;  // multi-species 2T

    // --- Outputs (owned by this class)
    surfaceVectorField phiUp_;      // full momentum flux
    surfaceScalarField phiEp_;      // total energy flux
    surfaceScalarField max_a_;      // max propagation speed (used in deltaT)

    // Vibrational outputs
    surfaceScalarField phiEv_;                 // convective flux of rho*eV
    PtrList<surfaceScalarField> phiRhoEvibYi_; // convective fluxes of rho*Yi*eVib_i

private:


public:
// ---------------------------------------------------------------------
    // --- 1) Directed interpolation of primitive fields onto faces
    // ---------------------------------------------------------------------
    void calcDirectedInterpolation()
    {
        rho_pos_ = interpolate(rho_, pos_);
        rho_neg_ = interpolate(rho_, neg_);

        rhoU_pos_ = interpolate(rhoU_, pos_, U_.name());
        rhoU_neg_ = interpolate(rhoU_, neg_, U_.name());

        // rPsi = 1/psi (EOS-dependent)
        volScalarField rPsi("rPsi", 1.0/psi_);
        rPsi_pos_ = interpolate(rPsi, pos_, TTR_.name());
        rPsi_neg_ = interpolate(rPsi, neg_, TTR_.name());

        e_pos_ = interpolate(eTR_, pos_, TTR_.name());
        e_neg_ = interpolate(eTR_, neg_, TTR_.name());

        U_pos_ = rhoU_pos_/rho_pos_;
        U_neg_ = rhoU_neg_/rho_neg_;

        p_pos_ = rho_pos_*rPsi_pos_;
        p_neg_ = rho_neg_*rPsi_neg_;

        // Face-normal volumetric fluxes (unoriented for max/min ops)
        phi_fp_ = (U_pos_ & mesh_.Sf());
        phi_fp_.setOriented(false);

        phi_fn_ = (U_neg_ & mesh_.Sf());
        phi_fn_.setOriented(false);

        // Make fluxes relative to mesh-motion
        if (mesh_.moving())
        {
            surfaceScalarField meshPhi(mesh_.phi());
            meshPhi.setOriented(false);
            phi_fp_ -= meshPhi;
            phi_fn_ -= meshPhi;
        }
    }

    // ---------------------------------------------------------------------
    // --- 2) Sound speed and face sound-speed*area
    // ---------------------------------------------------------------------
    void calcSoundSpeed()
    {
        // c = sqrt(gamma/psi)  with gamma ~ Cp/Cv (here TR mode)
        volScalarField rPsi("rPsi", 1.0/psi_);
        volScalarField c
        (
            "c",
            sqrt(thermo2T_.gamma()*rPsi)
        );

        cSf_pos_ = interpolate(c, pos_, TTR_.name())*mesh_.magSf();
        cSf_neg_ = interpolate(c, neg_, TTR_.name())*mesh_.magSf();
    }

    // ---------------------------------------------------------------------
    // --- 3) Local propagation speeds and KNP/KT weights
    // ---------------------------------------------------------------------
    void calcPropagationSpeedsAndWeights()
    {
        a_fp_ = max(max(phi_fp_ + cSf_pos_, phi_fn_ + cSf_neg_), v_zero_);
        a_fn_ = min(min(phi_fp_ - cSf_pos_, phi_fn_ - cSf_neg_), v_zero_);

        // Weighting factor, default KNP
        alpha_ = a_fp_/(a_fp_ - a_fn_);
        minus_alpha_ = 1.0 - alpha_;

        // max(|a-|,|a+|) used by KT/Tadmor and also for dt estimates
        max_a_ = max(mag(a_fn_), mag(a_fp_));

        // Diffusive flux factor, default KNP
        nw_f_ = a_fn_*alpha_;

        // KT/Tadmor corrections
        if (fluxScheme_ == "Tadmor")
        {
            alpha_ = 0.5;
            minus_alpha_ = 0.5;
            nw_f_ = -0.5*max_a_;
        }
    }

    // ---------------------------------------------------------------------
    // --- 4) Build final aphi coefficients and update max_a
    // ---------------------------------------------------------------------
    void calcAphi()
    {
        // Weighted flux calculations
        // phi_f+ * alpha , phi_f- * (1 - alpha)
        surfaceScalarField phi_fp_w(phi_fp_*alpha_);
        surfaceScalarField phi_fn_w(phi_fn_*minus_alpha_);

        // Diffusive corrections
        aphi_fp_ = phi_fp_w - nw_f_;
        aphi_fn_ = phi_fn_w + nw_f_;

        // Maximum flux magnitude used in deltaT calculations
        max_a_ = max(mag(aphi_fp_), mag(aphi_fn_));
    }

    // ---------------------------------------------------------------------
    // --- 5a) Face fluxes for conservative variables
    // ---------------------------------------------------------------------
    // -----------------------------------------------------------------------
    // Mass flux
    // -----------------------------------------------------------------------
    void calcMassFlux()
    {
        phi_ = aphi_fp_*rho_pos_ + aphi_fn_*rho_neg_;
        phi_.setOriented(true);
    }

    // -----------------------------------------------------------------------
    // Momentum flux (convective + pressure)
    // -----------------------------------------------------------------------
    void calcMomentumFlux()
    {
        surfaceVectorField phiU(aphi_fp_*rhoU_pos_ + aphi_fn_*rhoU_neg_);
        phiU.setOriented(true);

        phiUp_ =
            phiU
        + (alpha_*p_pos_ + minus_alpha_*p_neg_)*mesh_.Sf();
        phiUp_.setOriented(true);
    }

    // -----------------------------------------------------------------------
    // Total energy flux
    // -----------------------------------------------------------------------
    void calcEnergyFlux()
    {
        // --- Build face vibrational energy (zero if 1T) ---
        surfaceScalarField eVmix_pos
        (
            IOobject("eVmix_pos", mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar("zero", eTR_.dimensions(), 0.0)
        );

        surfaceScalarField eVmix_neg
        (
            IOobject("eVmix_neg", mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar("zero", eTR_.dimensions(), 0.0)
        );

        if (use2T_)
        {
            if (!mixtureCheck_ && eVPtr_)
            {
                eVmix_pos = interpolate(*eVPtr_, pos_, TTR_.name());
                eVmix_neg = interpolate(*eVPtr_, neg_, TTR_.name());
            }
            else if (mixtureCheck_ && eVibSpeciesPtr_)
            {
                const PtrList<volScalarField>& Y      = thermo2T_.composition().Y();
                const PtrList<volScalarField>& eVibSp = *eVibSpeciesPtr_;

                forAll(Y, i)
                {
                    if (!thermo2T_.composition().isSpecieMolecular(i)) continue;

                    surfaceScalarField Yi_pos(interpolate(Y[i], pos_, TTR_.name()));
                    surfaceScalarField Yi_neg(interpolate(Y[i], neg_, TTR_.name()));

                    Yi_pos = max(min(Yi_pos, scalar(1)), scalar(0));
                    Yi_neg = max(min(Yi_neg, scalar(1)), scalar(0));

                    eVmix_pos += Yi_pos*interpolate(eVibSp[i], pos_, TTR_.name());
                    eVmix_neg += Yi_neg*interpolate(eVibSp[i], neg_, TTR_.name());
                }
            }
        }

        // --- Kurganov energy flux ---
        phiEp_ =
            aphi_fp_*(rho_pos_*(e_pos_ + eVmix_pos + 0.5*magSqr(U_pos_)) + p_pos_)
        + aphi_fn_*(rho_neg_*(e_neg_ + eVmix_neg + 0.5*magSqr(U_neg_)) + p_neg_)
        + nw_f_*(p_pos_ - p_neg_);

        phiEp_.setOriented(true);

        // --- ALE correction for moving meshes ---
        if (mesh_.moving())
        {
            surfaceScalarField meshPhi(mesh_.phi());
            meshPhi.setOriented(false);
            phiEp_ += meshPhi*(alpha_*p_pos_ + minus_alpha_*p_neg_);
        }
    }

    void calcConservativeFluxes()
    {
        calcMassFlux();
        calcMomentumFlux();
        calcEnergyFlux();
    }
    
    // ---------------------------------------------------------------------
    // --- 5b) Face fluxes for multispecies
    // ---------------------------------------------------------------------
    // -----------------------------------------------------------------------
    // Mass fraction fluxes
    // -----------------------------------------------------------------------
    void calcMassFractionFluxes()
    {
        if (!mixtureCheck_ || nSpecie_ == 0) return;

        const PtrList<volScalarField>& Y = thermo2T_.composition().Y();
        const label nSpec = Y.size();
        const label last  = nSpec - 1;

        surfaceScalarField sumPhiRhoYi(phi_);
        sumPhiRhoYi.rename("sumPhiRhoYi");
        sumPhiRhoYi = dimensionedScalar(sumPhiRhoYi.dimensions(), 0.0);

        for (label i = 0; i < last; ++i)
        {
            surfaceScalarField Yi_pos(interpolate(Y[i], pos_, TTR_.name()));
            surfaceScalarField Yi_neg(interpolate(Y[i], neg_, TTR_.name()));

            Yi_pos = max(min(Yi_pos, scalar(1)), scalar(0));
            Yi_neg = max(min(Yi_neg, scalar(1)), scalar(0));

            const surfaceScalarField rhoYi_pos(rho_pos_*Yi_pos);
            const surfaceScalarField rhoYi_neg(rho_neg_*Yi_neg);

            phiRhoYi_[i] = aphi_fp_*rhoYi_pos + aphi_fn_*rhoYi_neg;
            phiRhoYi_[i].setOriented(true);

            sumPhiRhoYi += phiRhoYi_[i];
        }

        // Last species flux — closure to enforce sum(Y) = 1 exactly
        phiRhoYi_[last] = phi_ - sumPhiRhoYi;
        phiRhoYi_[last].setOriented(true);
    }

    void calcSpeciesFluxes()
    {
        calcMassFlux();
        calcMomentumFlux();
        calcEnergyFlux();       
        calcMassFractionFluxes();
    }

    void calcPureMixtureVibFlux()
    {
        if (!use2T_ || mixtureCheck_ || !eVPtr_)
        {
            return;
        }

        const volScalarField& eV = *eVPtr_;

        surfaceScalarField eV_pos(interpolate(eV, pos_, TTR_.name()));
        surfaceScalarField eV_neg(interpolate(eV, neg_, TTR_.name()));

        phiEv_ =
            aphi_fp_*(rho_pos_*eV_pos)
        + aphi_fn_*(rho_neg_*eV_neg);

        phiEv_.setOriented(true);
    }

    void calcSpeciesVibFluxes()
    {
        if (!use2T_ || !mixtureCheck_ || !eVibSpeciesPtr_ || nSpecie_ == 0)
        {
            return;
        }

        const PtrList<volScalarField>& Y = thermo2T_.composition().Y();
        const PtrList<volScalarField>& eVibSpecies = *eVibSpeciesPtr_;

        forAll(Y, i)
        {
            if (!thermo2T_.composition().isSpecieMolecular(i))
            {
                phiRhoEvibYi_[i] =
                    dimensionedScalar("zero", dimEnergy/dimTime, 0.0);
                phiRhoEvibYi_[i].setOriented(true);
                continue;
            }

            surfaceScalarField Yi_pos(interpolate(Y[i], pos_, TTR_.name()));
            surfaceScalarField Yi_neg(interpolate(Y[i], neg_, TTR_.name()));

            Yi_pos = max(min(Yi_pos, scalar(1)), scalar(0));
            Yi_neg = max(min(Yi_neg, scalar(1)), scalar(0));

            surfaceScalarField eV_pos
            (
                interpolate(eVibSpecies[i], pos_, TTR_.name())
            );

            surfaceScalarField eV_neg
            (
                interpolate(eVibSpecies[i], neg_, TTR_.name())
            );

            phiRhoEvibYi_[i] =
                aphi_fp_*(rho_pos_*Yi_pos*eV_pos)
            + aphi_fn_*(rho_neg_*Yi_neg*eV_neg);

            phiRhoEvibYi_[i].setOriented(true);
        }
    }
    
    /***************************** Constructors ***************************/
    // Existing/base path: 1T or "no vibrational flux requested"
    kurganovFluxCalculator
    (
        const fvMesh& mesh,
        tcLibraryInterface& thermo2T,
        const volScalarField& rho,
        const volVectorField& rhoU,
        const volVectorField& U,
        const volScalarField& psi,
        const volScalarField& eTR,
        const volScalarField& TTR,
        const surfaceScalarField& pos,
        const surfaceScalarField& neg,
        surfaceScalarField& phi,
        const word& fluxScheme,
        const dimensionedScalar& v_zero,
        const bool mixtureCheck,
        const bool use2T
    )
    :
        mesh_(mesh),
        thermo2T_(thermo2T),

        rho_(rho),
        rhoU_(rhoU),
        U_(U),

        psi_(psi),
        eTR_(eTR),
        TTR_(TTR),

        pos_(pos),
        neg_(neg),

        phi_(phi),

        fluxScheme_(fluxScheme),
        v_zero_(v_zero),

        mixtureCheck_(mixtureCheck),
        use2T_(use2T),
        
        eVPtr_(nullptr),
        eVibSpeciesPtr_(nullptr),

        // Workspace fields
        rho_pos_
        (
            IOobject
            (
                "rho_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions(), 
                0.0
            )
        ),

        rho_neg_
        (
            IOobject
            (
                "rho_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions(), 
                0.0
            )
        ),

        rhoU_pos_
        (
            IOobject
            (
                "rhoU_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions(), 
                vector::zero
            )
        ),

        rhoU_neg_
        (
            IOobject
            (
                "rhoU_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions(), 
                vector::zero
            )
        ),

        rPsi_pos_
        (
            IOobject
            (
                "rPsi_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure/dimDensity, 
                0.0
            )
        ),
        
        rPsi_neg_
        (
            IOobject
            (
                "rPsi_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure/dimDensity, 
                0.0
            )
        ),

        e_pos_
        (
            IOobject
            (
                "h_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                eTR_.dimensions(), 
                0.0
            )
        ),
        
        e_neg_
        (
            IOobject
            (
                "h_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                eTR_.dimensions(), 
                0.0
            )
        ),

        U_pos_
        (
            IOobject
            (
                "U_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                U_.dimensions(), 
                vector::zero
            )
        ),
        
        U_neg_
        (
            IOobject
            (
                "U_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                U_.dimensions(), 
                vector::zero
            )
        ),

        p_pos_
        (
            IOobject
            (
                "p_pos", 
                mesh_.time().timeName(), 
                mesh_, IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure, 
                0.0
            )
        ),
        
        p_neg_
        (
            IOobject
            (
                "p_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure, 
                0.0
            )
        ),

        phi_fp_
        (
            IOobject
            (
                "phi_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        phi_fn_
        (
            IOobject
            (
                "phi_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        cSf_pos_
        (
            IOobject
            (
                "cSf_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        cSf_neg_
        (
            IOobject
            (
                "cSf_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        a_fp_
        (
            IOobject
            (
                "a_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        a_fn_
        (
            IOobject
            (
                "a_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        alpha_
        (
            IOobject
            (
                "alpha", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimless, 
                0.0
            )
        ),

        minus_alpha_
        (
            IOobject
            (
                "minus_alpha", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimless, 
                0.0
            )
        ),

        nw_f_
        (
            IOobject
            (
                "nw_f", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        aphi_fp_
        (
            IOobject
            (
                "aphi_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        aphi_fn_
        (
            IOobject
            (
                "aphi_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        // Outputs
        phiUp_
        (
            IOobject
            (
                "phiUp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions()*dimVolume/dimTime, 
                vector::zero
            )
        ),

        phiEp_
        (
            IOobject
            (
                "phiEp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions()*eTR_.dimensions()*dimVolume/dimTime, 
                0.0
            )
        ),
        
        max_a_
        (
            IOobject
            (
                "max_a", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        phiEv_
        (
            IOobject
            (
                "phiEv",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
        )
    {
        if (mixtureCheck_)
        {
            // Get number of species from thermo
            // In reacting/mixture thermos, composition().Y() is the canonical list.
            // This compiles for standard OpenFOAM thermos that expose composition().
            const PtrList<volScalarField>& Y = thermo2T_.composition().Y();
            nSpecie_ = Y.size();

            phiRhoYi_.setSize(nSpecie_);

            forAll(phiRhoYi_, i)
            {
                // Name like phiRhoY_N2, phiRhoY_O2, ...
                const word fldName("phiRhoY_" + Y[i].name());

                phiRhoYi_.set
                (
                    i,
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            fldName,
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE  // or NO_WRITE if you don’t want files
                        ),
                        mesh_,
                        dimensionedScalar(dimMass/dimTime, 0.0) // kg/s
                    )
                );
                phiRhoYi_[i].setOriented(true);
            }
        }
        else
        {
            nSpecie_ = 0;
            phiRhoYi_.clear();
            phiRhoEvibYi_.clear();
        }
    }

    // 2T pure-mixture path: convective flux of rho*eV
    kurganovFluxCalculator
    (
        const fvMesh& mesh,
        tcLibraryInterface& thermo2T,
        const volScalarField& rho,
        const volVectorField& rhoU,
        const volVectorField& U,
        const volScalarField& psi,
        const volScalarField& eTR,
        const volScalarField& TTR,
        const volScalarField& eV,
        const surfaceScalarField& pos,
        const surfaceScalarField& neg,
        surfaceScalarField& phi,
        const word& fluxScheme,
        const dimensionedScalar& v_zero,
        const bool mixtureCheck,
        const bool use2T
    )
    :
        mesh_(mesh),
        thermo2T_(thermo2T),

        rho_(rho),
        rhoU_(rhoU),
        U_(U),

        psi_(psi),
        eTR_(eTR),
        TTR_(TTR),

        pos_(pos),
        neg_(neg),

        phi_(phi),

        fluxScheme_(fluxScheme),
        v_zero_(v_zero),

        mixtureCheck_(mixtureCheck),
        use2T_(use2T),
        
        eVPtr_(nullptr),
        eVibSpeciesPtr_(nullptr),

        // Workspace fields
        rho_pos_
        (
            IOobject
            (
                "rho_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions(), 
                0.0
            )
        ),

        rho_neg_
        (
            IOobject
            (
                "rho_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions(), 
                0.0
            )
        ),

        rhoU_pos_
        (
            IOobject
            (
                "rhoU_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions(), 
                vector::zero
            )
        ),

        rhoU_neg_
        (
            IOobject
            (
                "rhoU_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions(), 
                vector::zero
            )
        ),

        rPsi_pos_
        (
            IOobject
            (
                "rPsi_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure/dimDensity, 
                0.0
            )
        ),
        
        rPsi_neg_
        (
            IOobject
            (
                "rPsi_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure/dimDensity, 
                0.0
            )
        ),

        e_pos_
        (
            IOobject
            (
                "h_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                eTR_.dimensions(), 
                0.0
            )
        ),
        
        e_neg_
        (
            IOobject
            (
                "h_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                eTR_.dimensions(), 
                0.0
            )
        ),

        U_pos_
        (
            IOobject
            (
                "U_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                U_.dimensions(), 
                vector::zero
            )
        ),
        
        U_neg_
        (
            IOobject
            (
                "U_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                U_.dimensions(), 
                vector::zero
            )
        ),

        p_pos_
        (
            IOobject
            (
                "p_pos", 
                mesh_.time().timeName(), 
                mesh_, IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure, 
                0.0
            )
        ),
        
        p_neg_
        (
            IOobject
            (
                "p_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure, 
                0.0
            )
        ),

        phi_fp_
        (
            IOobject
            (
                "phi_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        phi_fn_
        (
            IOobject
            (
                "phi_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        cSf_pos_
        (
            IOobject
            (
                "cSf_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        cSf_neg_
        (
            IOobject
            (
                "cSf_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        a_fp_
        (
            IOobject
            (
                "a_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        a_fn_
        (
            IOobject
            (
                "a_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        alpha_
        (
            IOobject
            (
                "alpha", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimless, 
                0.0
            )
        ),

        minus_alpha_
        (
            IOobject
            (
                "minus_alpha", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimless, 
                0.0
            )
        ),

        nw_f_
        (
            IOobject
            (
                "nw_f", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        aphi_fp_
        (
            IOobject
            (
                "aphi_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        aphi_fn_
        (
            IOobject
            (
                "aphi_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        // Outputs
        phiUp_
        (
            IOobject
            (
                "phiUp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions()*dimVolume/dimTime, 
                vector::zero
            )
        ),

        phiEp_
        (
            IOobject
            (
                "phiEp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions()*eTR_.dimensions()*dimVolume/dimTime, 
                0.0
            )
        ),
        
        max_a_
        (
            IOobject
            (
                "max_a", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        phiEv_
        (
            IOobject
            (
                "phiEv",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
        )
    {
        eVPtr_ = &eV;
        phiEv_.setOriented(true);
        if (mixtureCheck_)
        {
            // Get number of species from thermo
            // In reacting/mixture thermos, composition().Y() is the canonical list.
            // This compiles for standard OpenFOAM thermos that expose composition().
            const PtrList<volScalarField>& Y = thermo2T_.composition().Y();
            nSpecie_ = Y.size();

            phiRhoYi_.setSize(nSpecie_);

            forAll(phiRhoYi_, i)
            {
                // Name like phiRhoY_N2, phiRhoY_O2, ...
                const word fldName("phiRhoY_" + Y[i].name());

                phiRhoYi_.set
                (
                    i,
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            fldName,
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE  // or NO_WRITE if you don’t want files
                        ),
                        mesh_,
                        dimensionedScalar(dimMass/dimTime, 0.0) // kg/s
                    )
                );
                phiRhoYi_[i].setOriented(true);
            }
        }
        else
        {
            nSpecie_ = 0;
            phiRhoYi_.clear();
            phiRhoEvibYi_.clear();
        }
    }

    // 2T multi-species path: convective fluxes of rho*Yi*eVib_i
    kurganovFluxCalculator
    (
        const fvMesh& mesh,
        tcLibraryInterface& thermo2T,
        const volScalarField& rho,
        const volVectorField& rhoU,
        const volVectorField& U,
        const volScalarField& psi,
        const volScalarField& eTR,
        const volScalarField& TTR,
        const PtrList<volScalarField>& eVibSpecies,
        const surfaceScalarField& pos,
        const surfaceScalarField& neg,
        surfaceScalarField& phi,
        const word& fluxScheme,
        const dimensionedScalar& v_zero,
        const bool mixtureCheck,
        const bool use2T
    )
    :
        mesh_(mesh),
        thermo2T_(thermo2T),

        rho_(rho),
        rhoU_(rhoU),
        U_(U),

        psi_(psi),
        eTR_(eTR),
        TTR_(TTR),

        pos_(pos),
        neg_(neg),

        phi_(phi),

        fluxScheme_(fluxScheme),
        v_zero_(v_zero),

        mixtureCheck_(mixtureCheck),
        use2T_(use2T),
        
        eVPtr_(nullptr),
        eVibSpeciesPtr_(nullptr),

        // Workspace fields
        rho_pos_
        (
            IOobject
            (
                "rho_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions(), 
                0.0
            )
        ),

        rho_neg_
        (
            IOobject
            (
                "rho_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions(), 
                0.0
            )
        ),

        rhoU_pos_
        (
            IOobject
            (
                "rhoU_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions(), 
                vector::zero
            )
        ),

        rhoU_neg_
        (
            IOobject
            (
                "rhoU_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions(), 
                vector::zero
            )
        ),

        rPsi_pos_
        (
            IOobject
            (
                "rPsi_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure/dimDensity, 
                0.0
            )
        ),
        
        rPsi_neg_
        (
            IOobject
            (
                "rPsi_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure/dimDensity, 
                0.0
            )
        ),

        e_pos_
        (
            IOobject
            (
                "h_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                eTR_.dimensions(), 
                0.0
            )
        ),
        
        e_neg_
        (
            IOobject
            (
                "h_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                eTR_.dimensions(), 
                0.0
            )
        ),

        U_pos_
        (
            IOobject
            (
                "U_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                U_.dimensions(), 
                vector::zero
            )
        ),
        
        U_neg_
        (
            IOobject
            (
                "U_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                U_.dimensions(), 
                vector::zero
            )
        ),

        p_pos_
        (
            IOobject
            (
                "p_pos", 
                mesh_.time().timeName(), 
                mesh_, IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure, 
                0.0
            )
        ),
        
        p_neg_
        (
            IOobject
            (
                "p_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimPressure, 
                0.0
            )
        ),

        phi_fp_
        (
            IOobject
            (
                "phi_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        phi_fn_
        (
            IOobject
            (
                "phi_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        cSf_pos_
        (
            IOobject
            (
                "cSf_pos", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        cSf_neg_
        (
            IOobject
            (
                "cSf_neg", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        a_fp_
        (
            IOobject
            (
                "a_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        a_fn_
        (
            IOobject
            (
                "a_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        alpha_
        (
            IOobject
            (
                "alpha", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimless, 
                0.0
            )
        ),

        minus_alpha_
        (
            IOobject
            (
                "minus_alpha", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimless, 
                0.0
            )
        ),

        nw_f_
        (
            IOobject
            (
                "nw_f", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        aphi_fp_
        (
            IOobject
            (
                "aphi_fp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),
        
        aphi_fn_
        (
            IOobject
            (
                "aphi_fn", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        // Outputs
        phiUp_
        (
            IOobject
            (
                "phiUp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedVector
            (
                "zero", 
                rhoU_.dimensions()*dimVolume/dimTime, 
                vector::zero
            )
        ),

        phiEp_
        (
            IOobject
            (
                "phiEp", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                rho_.dimensions()*eTR_.dimensions()*dimVolume/dimTime, 
                0.0
            )
        ),
        
        max_a_
        (
            IOobject
            (
                "max_a", 
                mesh_.time().timeName(), 
                mesh_, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            mesh_, 
            dimensionedScalar
            (
                "zero", 
                dimVolume/dimTime, 
                0.0
            )
        ),

        phiEv_
        (
            IOobject
            (
                "phiEv",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
        )
        {
        eVibSpeciesPtr_ = &eVibSpecies;
        phiEv_.setOriented(true);

        if (mixtureCheck_)
        {
            const PtrList<volScalarField>& Y = thermo2T_.composition().Y();
            nSpecie_ = Y.size();

            phiRhoYi_.setSize(nSpecie_);
            phiRhoEvibYi_.setSize(nSpecie_);

            forAll(Y, i)
            {
                const word fldNameMass("phiRhoY_" + Y[i].name());
                const word fldNameVib("phiRhoEvib_" + Y[i].name());

                phiRhoYi_.set
                (
                    i,
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            fldNameMass,
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar(dimMass/dimTime, 0.0)
                    )
                );
                phiRhoYi_[i].setOriented(true);

                phiRhoEvibYi_.set
                (
                    i,
                    new surfaceScalarField
                    (
                        IOobject
                        (
                            fldNameVib,
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar(dimEnergy/dimTime, 0.0)
                    )
                );
                phiRhoEvibYi_[i].setOriented(true);
            }
        }
        else
        {
            nSpecie_ = 0;
            phiRhoYi_.clear();
            phiRhoEvibYi_.clear();
        }
    }

    // Accessors for solver usage
    inline surfaceVectorField&          phiUp()          { return phiUp_; }
    inline surfaceScalarField&          phiEp()          { return phiEp_; }
    inline const  PtrList<surfaceScalarField>& phiRhoYi() const { return phiRhoYi_; }
    inline surfaceScalarField&          maxA()           { return max_a_;  }
    inline surfaceScalarField& phiEv() { return phiEv_; }
    inline const PtrList<surfaceScalarField>& phiRhoEvibYi() const
    {
        return phiRhoEvibYi_;
    }
    inline surfaceVectorField&          U_pos()          { return U_pos_; }
    inline surfaceVectorField&          U_neg()          { return U_neg_; }
    // Optional: access alpha/nw_f for diagnostics
    inline surfaceScalarField& alpha() { return alpha_; }
    inline surfaceScalarField& minus_alpha() { return minus_alpha_; }
    inline const surfaceScalarField& nw_f()  const { return nw_f_;  }
};

} // namespace Foam

#endif
