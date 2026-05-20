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

#include "smoluchowskiJumpT.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "tcLibraryInterface.H" 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoluchowskiJumpT::smoluchowskiJumpT
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho"),
    psiName_("thermo:psi"),
    muName_("thermo:mu"),
    accommodationCoeff_(1.0),
    Twall_(p.size(), Zero)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::smoluchowskiJumpT::smoluchowskiJumpT
(
    const smoluchowskiJumpT& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    muName_(ptf.muName_),
    accommodationCoeff_(ptf.accommodationCoeff_),
    Twall_(ptf.Twall_)
{}


Foam::smoluchowskiJumpT::smoluchowskiJumpT
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.getOrDefault<word>("U", "U")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "thermo:psi")),
    muName_(dict.getOrDefault<word>("mu", "thermo:mu")),
    accommodationCoeff_(dict.get<scalar>("accommodationCoeff")),
    Twall_("Twall", dict, p.size())
{
    if
    (
        mag(accommodationCoeff_) < SMALL
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorInFunction(dict)
            << "unphysical accommodationCoeff specified"
            << "(0 < accommodationCoeff <= 1)" << endl
            << exit(FatalIOError);
    }

    if (!this->readValueEntry(dict))
    {
        // Fallback: set to the internal field
        fvPatchField<scalar>::patchInternalField(*this);
    }

    refValue() = *this;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::smoluchowskiJumpT::smoluchowskiJumpT
(
    const smoluchowskiJumpT& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    UName_(ptpsf.UName_),
    rhoName_(ptpsf.rhoName_),
    psiName_(ptpsf.psiName_),
    muName_(ptpsf.muName_),
    accommodationCoeff_(ptpsf.accommodationCoeff_),
    Twall_(ptpsf.Twall_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoluchowskiJumpT::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    Twall_.autoMap(m);
}


void Foam::smoluchowskiJumpT::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const smoluchowskiJumpT& tiptf =
        refCast<const smoluchowskiJumpT>(ptf);

    Twall_.rmap(tiptf.Twall_, addr);
}


void Foam::smoluchowskiJumpT::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Guard: skip if thermo stack is not yet fully initialised
    if
    (
       !db().foundObject<volScalarField>("thermo:alpha")
    || !db().foundObject<volScalarField>(muName_)
    || !db().foundObject<volScalarField>(psiName_)
    || !db().foundObject<volScalarField>(rhoName_)
    || !db().foundObject<tcLibraryInterface>("thermophysicalProperties")
    )
    {
        valueFraction() = 0.0;
        refValue() = Twall_;
        refGrad() = 0.0;
        return;
    }

    // --- Field lookups from thermo stack
    const auto& pmu  = patch().lookupPatchField<volScalarField>(muName_);
    const auto& prho = patch().lookupPatchField<volScalarField>(rhoName_);
    const auto& ppsi = patch().lookupPatchField<volScalarField>(psiName_);

    // --- Thermal diffusivity alpha = kappaTR/Cp — registered by thermo stack
    const auto& palpha =
        patch().lookupPatchField<volScalarField>("thermo:alpha");

    // --- Prandtl number: face-varying, derived from thermo stack
    //
    //     Pr = mu / alpha_TR
    //
    //     thermo:alpha = kappaTR/Cp is already registered — no separate
    //     kappaTR or Cp lookup needed. Exact for any mixture.
    const Field<scalar> Pr(pmu/palpha);

    // --- Ratio of specific heats: face-varying, from tcLibraryInterface
    //
    //     gamma() is a virtual method on tcLibraryInterface — exact,
    //     computed by the thermo backend for the current composition
    //     and temperature. No fixed scalar assumption.
    const tcLibraryInterface& thermo =
        db().lookupObject<tcLibraryInterface>("thermophysicalProperties");

    const label patchi = patch().index();
    const Field<scalar> gamma
    (
        thermo.gamma()().boundaryField()[patchi]
    );

    // --- Smoluchowski jump length scale C2 (Eq. 127 prefactor)
    //
    //     C2 = (mu/rho) * sqrt(psi*pi/2) * 2*gamma/((gamma+1)*Pr)
    //          * (2 - alpha) / alpha
    //
    //     All quantities are face-varying. No calorically perfect assumption.
    Field<scalar> C2
    (
        pmu/prho
      * sqrt(ppsi*constant::mathematical::piByTwo)
      * 2.0*gamma/Pr/(gamma + 1.0)
      * (2.0 - accommodationCoeff_)/accommodationCoeff_
    );

    // --- Robin blending coefficient
    //
    //     valueFraction -> 1: Kn -> 0, hTR_patch -> hTRwall (no-jump)
    //     valueFraction -> 0: Kn -> inf, TTR drifts to interior value
    valueFraction() = 1.0/(1.0 + patch().deltaCoeffs()*C2);

    // --- Reference value: wall temperature directly
    //
    //     fixesValue() = true causes calculate() to back-compute
    //     hTR from this TTR value — exact, no Cp multiplication needed.
    refValue() = Twall_;

    refGrad() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::smoluchowskiJumpT::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    os.writeEntryIfDifferent<word>("mu", "thermo:mu", muName_);

    os.writeEntry("accommodationCoeff", accommodationCoeff_);
    Twall_.writeEntry("Twall", os);

    // Note: gamma and Pr are not written — they are derived from the
    // thermo stack at runtime and must not be stored as fixed scalars.

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        smoluchowskiJumpT
    );
}

// ************************************************************************* //