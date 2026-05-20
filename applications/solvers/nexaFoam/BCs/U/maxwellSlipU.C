/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "maxwellSlipU.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::maxwellSlipU::maxwellSlipU
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    TName_("TTR"),
    rhoName_("rho"),
    psiName_("thermo2T:psi"),
    muName_("thermo2T:mu"),
    tauMCName_("tauMC"),
    accommodationCoeff_(1.0),
    Uwall_(p.size(), Zero),
    thermalCreep_(true),
    curvature_(true)
{}


Foam::maxwellSlipU::maxwellSlipU
(
    const maxwellSlipU& mspvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(mspvf, p, iF, mapper),
    TName_(mspvf.TName_),
    rhoName_(mspvf.rhoName_),
    psiName_(mspvf.psiName_),
    muName_(mspvf.muName_),
    tauMCName_(mspvf.tauMCName_),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


Foam::maxwellSlipU::maxwellSlipU
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    TName_(dict.getOrDefault<word>("T", "TTR")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "thermo2T:psi")),
    muName_(dict.getOrDefault<word>("mu", "thermo2T:mu")),
    tauMCName_(dict.getOrDefault<word>("tauMC", "tauMC")),
    accommodationCoeff_(dict.get<scalar>("accommodationCoeff")),
    Uwall_("Uwall", dict, p.size()),
    thermalCreep_(dict.getOrDefault("thermalCreep", true)),
    curvature_(dict.getOrDefault("curvature", true))
{
    if
    (
        mag(accommodationCoeff_) < SMALL
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorInFunction(dict)
            << "unphysical accommodationCoeff_ specified"
            << "(0 < accommodationCoeff_ <= 1)" << endl
            << exit(FatalIOError);
    }

    if (this->readValueEntry(dict))
    {
        const auto* hasRefValue = dict.findEntry("refValue", keyType::LITERAL);
        const auto* hasFrac     = dict.findEntry("valueFraction", keyType::LITERAL);

        if (hasRefValue && hasFrac)
        {
            this->refValue().assign(*hasRefValue, p.size());
            this->valueFraction().assign(*hasFrac, p.size());
        }
        else
        {
            this->refValue() = *this;
            this->valueFraction() = scalar(1);
        }
    }
}


Foam::maxwellSlipU::maxwellSlipU
(
    const maxwellSlipU& mspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(mspvf, iF),
    TName_(mspvf.TName_),
    rhoName_(mspvf.rhoName_),
    psiName_(mspvf.psiName_),
    muName_(mspvf.muName_),
    tauMCName_(mspvf.tauMCName_),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::maxwellSlipU::autoMap
(
    const fvPatchFieldMapper& m
)
{
    partialSlipFvPatchVectorField::autoMap(m);
    Uwall_.autoMap(m);
}


void Foam::maxwellSlipU::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    partialSlipFvPatchVectorField::rmap(ptf, addr);

    const maxwellSlipU& mspvf =
        refCast<const maxwellSlipU>(ptf);

    Uwall_.rmap(mspvf.Uwall_, addr);
}


void Foam::maxwellSlipU::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if
    (
       !db().foundObject<volScalarField>(muName_)
    || !db().foundObject<volScalarField>(rhoName_)
    || !db().foundObject<volScalarField>(psiName_)
    )
    {
        valueFraction() = scalar(1);
        refValue() = Uwall_;
        return;
    }

    // --- Field lookups from thermo stack
    const auto& pmu  = patch().lookupPatchField<volScalarField>(muName_);
    const auto& prho = patch().lookupPatchField<volScalarField>(rhoName_);
    const auto& ppsi = patch().lookupPatchField<volScalarField>(psiName_);

    // --- Wall normal vector
    //
    //     Factored out here rather than recomputed inside each optional block.
    //     Used by both the thermal creep and curvature corrections to project
    //     vector quantities onto the wall-tangent plane via transform(I - n*n, .)
    const vectorField n(patch().nf());

    // --- Kinematic viscosity
    //
    //     nu = mu / rho, derived from the thermo stack at the wall faces.
    //     Used in both the slip fraction and the thermal creep prefactor.
    const Field<scalar> pnu(pmu/prho);

    // --- Maxwell slip length scale C1
    //
    //     C1 = sqrt(psi*pi/2) * (2 - alpha) / alpha
    //
    //     where:
    //       sqrt(psi*pi/2)  =  lambda/nu  (mean free path / kinematic viscosity,
    //                          from kinetic theory Eq. 128 with psi = 1/(R*T))
    //       (2-alpha)/alpha =  Maxwell momentum accommodation factor (Eq. 126)
    //
    //     C1 is face-varying because psi varies along the wall with local
    //     temperature and composition.
    const Field<scalar> C1
    (
        sqrt(ppsi*constant::mathematical::piByTwo)
      * (2.0 - accommodationCoeff_)/accommodationCoeff_
    );

    // --- Robin blending coefficient
    //
    //     valueFraction = 1 / (1 + delta * C1 * nu)
    //
    //     Controls the weighting between refValue (wall velocity + corrections)
    //     and the interior velocity. Limits:
    //       Kn -> 0  (C1 -> 0): valueFraction -> 1, U_patch -> Uwall (no-slip)
    //       Kn -> inf (C1 large): valueFraction -> 0, U drifts to interior
    //
    //     The partialSlip base class additionally enforces zero wall-normal
    //     velocity regardless of valueFraction.
    valueFraction() = 1.0/(1.0 + patch().deltaCoeffs()*C1*pnu);

    // --- Base reference value: wall velocity
    //
    //     Thermal creep and curvature corrections are added on top of this.
    //     For a stationary wall Uwall_ = (0 0 0) everywhere.
    refValue() = Uwall_;

    // --- Thermal creep correction (optional, default on)
    //
    //     refValue -= 3*nu/(4*TTR) * tangentialComponent(grad(TTR))
    //
    //     Implements the sigma_T term from Maxwell's condition (Eq. 125).
    //     A tangential temperature gradient along the wall drives a slip
    //     velocity from cold toward hot regions independently of shear.
    //     Only the wall-parallel component of grad(TTR) contributes —
    //     transform(I - n*n, .) removes the normal component.
    //
    //     TTR is looked up by TName_ which defaults to "TTR" — the
    //     translational-rotational temperature field in nexaFoam. T_vib
    //     does not enter this term; creep is driven by translational
    //     kinetic energy gradients only.
    //
    //     Switch off with thermalCreep false; in the BC dictionary.
    if (thermalCreep_)
    {
        const volScalarField& vsfT =
            this->db().objectRegistry::lookupObject<volScalarField>(TName_);

        const label patchi = this->patch().index();
        const fvPatchScalarField& pT = vsfT.boundaryField()[patchi];
        const Field<vector> gradpT(fvc::grad(vsfT)().boundaryField()[patchi]);

        refValue() -= 3.0*pnu/(4.0*pT)*transform(I - n*n, gradpT);
    }

    // --- Surface curvature correction (optional, default on)
    //
    //     refValue -= C1/rho * tangentialComponent(n . tauMC)
    //
    //     Accounts for the tangential traction from the non-equilibrium
    //     Maxwell stress tensor on curved surfaces. tauMC is the deviatoric
    //     viscous stress tensor computed in nexaFoam.C as:
    //
    //         tauMC = muLam * dev2(T(fvc::grad(U)))
    //
    //     (n & tauMC) contracts the tensor with the wall normal to extract
    //     the traction vector. transform(I - n*n, .) projects it onto the
    //     tangent plane. tauMC must be registered in the object registry
    //     before U.correctBoundaryConditions() is called each time step.
    //
    //     Switch off with curvature false; in the BC dictionary.
    if (curvature_)
    {
        const auto& ptauMC =
            patch().lookupPatchField<volTensorField>(tauMCName_);

        refValue() -= C1/prho*transform(I - n*n, (n & ptauMC));
    }

    partialSlipFvPatchVectorField::updateCoeffs();
}


void Foam::maxwellSlipU::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeEntryIfDifferent<word>("T",     "TTR",          TName_);
    os.writeEntryIfDifferent<word>("rho",   "rho",          rhoName_);
    os.writeEntryIfDifferent<word>("psi",   "thermo2T:psi", psiName_);
    os.writeEntryIfDifferent<word>("mu",    "thermo2T:mu",  muName_);
    os.writeEntryIfDifferent<word>("tauMC", "tauMC",        tauMCName_);

    os.writeEntry("accommodationCoeff", accommodationCoeff_);
    Uwall_.writeEntry("Uwall", os);
    os.writeEntry("thermalCreep", thermalCreep_);
    os.writeEntry("curvature", curvature_);

    refValue().writeEntry("refValue", os);
    valueFraction().writeEntry("valueFraction", os);

    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        maxwellSlipU
    );
}

// ************************************************************************* //