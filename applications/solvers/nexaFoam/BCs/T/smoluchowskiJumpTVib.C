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

#include "smoluchowskiJumpTVib.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "tcLibraryInterface.H" 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoluchowskiJumpTVib::smoluchowskiJumpTVib
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    UName_("U"),
    rhoName_("rho"),
    psiName_("thermo2T:psi"),
    muName_("thermo2T:mu"),
    accommodationCoeffVib_(1.0),
    TVibWall_(p.size(), Zero)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::smoluchowskiJumpTVib::smoluchowskiJumpTVib
(
    const smoluchowskiJumpTVib& ptf,
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
    accommodationCoeffVib_(ptf.accommodationCoeffVib_),
    TVibWall_(ptf.TVibWall_)
{}


Foam::smoluchowskiJumpTVib::smoluchowskiJumpTVib
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    UName_(dict.getOrDefault<word>("U", "U")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "thermo2T:psi")),
    muName_(dict.getOrDefault<word>("mu", "thermo2T:mu")),
    accommodationCoeffVib_(dict.get<scalar>("accommodationCoeffVib")),
    TVibWall_("TVibWall", dict, p.size())
{
    if
    (
        mag(accommodationCoeffVib_) < SMALL
     || mag(accommodationCoeffVib_) > 2.0
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


Foam::smoluchowskiJumpTVib::smoluchowskiJumpTVib
(
    const smoluchowskiJumpTVib& ptpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptpsf, iF),
    UName_(ptpsf.UName_),
    rhoName_(ptpsf.rhoName_),
    psiName_(ptpsf.psiName_),
    muName_(ptpsf.muName_),
    accommodationCoeffVib_(ptpsf.accommodationCoeffVib_),
    TVibWall_(ptpsf.TVibWall_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoluchowskiJumpTVib::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    TVibWall_.autoMap(m);
}


void Foam::smoluchowskiJumpTVib::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const smoluchowskiJumpTVib& tiptf =
        refCast<const smoluchowskiJumpTVib>(ptf);

    TVibWall_.rmap(tiptf.TVibWall_, addr);
}


void Foam::smoluchowskiJumpTVib::updateCoeffs()
{
    if (updated()) return;

    if (!db().foundObject<tcLibraryInterface>("thermophysicalProperties"))
    {
        valueFraction() = 0.0;
        refValue() = TVibWall_;
        refGrad() = 0.0;
        return;
    }

    const tcLibraryInterface& thermo =
        db().lookupObject<tcLibraryInterface>("thermophysicalProperties");

    const label patchi = patch().index();
    const labelUList& fc = patch().faceCells();

    const tmp<volScalarField> tmu  = thermo.mu();
    const scalarField& muInt       = tmu()();
    const scalarField& psiInt      = thermo.psi()();
    const tmp<volScalarField> trho = thermo.rho();
    const scalarField& rhoInt      = trho()();
    const scalarField& TVibInt =
        db().lookupObject<volScalarField>("TVib")();     // mixture vib temperature (cells)

    // --- vibrational diffusivity alphaVe = kappaVib / CvVib  (Pr_vib = mu/alphaVe) ---
    // DECISION POINT: pick ONE of:
    //   (A) if the thermo exposes it directly:
    //         const scalarField& alphaVeInt = thermo.alphaheVib()();
    //   (B) else build it from kappaVib and CvVib (as vibEnergyEquation.H does):
    //         alphaVe = thermo.kappaVib() / CvVibMix    (per near-wall cell)
    const tmp<volScalarField> tkVib = thermo.kappaVib();
    const scalarField& kVibInt = tkVib()();

    Field<scalar> pmu(patch().size()), ppsi(patch().size()),
                  prho(patch().size()), pkVib(patch().size());
    forAll(fc, i)
    {
        pmu[i]   = muInt[fc[i]];
        ppsi[i]  = psiInt[fc[i]];
        prho[i]  = rhoInt[fc[i]];
        pkVib[i] = kVibInt[fc[i]];
    }

    const Field<scalar> gamma(thermo.gamma()().boundaryField()[patchi]);

    // --- Vibrational Prandtl number: Pr_vib = mu * Cv_vib / kappa_vib -------
    // Cv_vib mirrors vibEnergyEquation.H CvVibMix (SHO, molecular species only,
    // mass-fraction weighted) so the jump BC and the vib-conduction solve
    // use an IDENTICAL vibrational specific heat.
    forAll(fc, i) pkVib[i] = kVibInt[fc[i]];

    // Species data needed for the SHO Cv_vib
    const compositionInterface& comp = thermo.composition();
    const label nSpec = comp.Y().size();

    // Cv_vib per near-wall cell (mass-specific, J/kg/K)
    Field<scalar> pCvVib(patch().size(), Zero);
    forAll(fc, i)
    {
        const label celli = fc[i];
        scalar sumc = 0.0;
        for (label s = 0; s < nSpec; ++s)
        {
            if (!comp.isSpecieMolecular(s)) continue;
            const scalar thetaV = comp.thetaVib(s);
            if (thetaV < SMALL) continue;

            const scalar Tv = max(TVibInt[celli], scalar(1));   // ← was comp.TVibSpecies(s)
            const scalar x  = 0.5*thetaV/Tv;
            const scalar xs = (x < 20.0) ? x/Foam::sinh(x) : 0.0;

            sumc += comp.Y()[s][celli] * comp.R(s) * xs*xs;
        }
        pCvVib[i] = max(sumc, SMALL);
    }

    // Pr_vib, guarded against zero kappaVib
    Field<scalar> Pr_vib(patch().size());
    forAll(Pr_vib, i)
    {
        Pr_vib[i] = (mag(pkVib[i]) > SMALL)
                  ? pmu[i]*pCvVib[i]/pkVib[i]
                  : 1.0;                                       // safe fallback
    }

    Field<scalar> C2
    (
        pmu/prho
      * sqrt(ppsi*constant::mathematical::piByTwo)
      * 2.0*gamma/Pr_vib/(gamma + 1.0)
      * (2.0 - accommodationCoeffVib_)/accommodationCoeffVib_
    );

    valueFraction() = 1.0/(1.0 + patch().deltaCoeffs()*C2);
    refValue()      = TVibWall_;   // full accommodation: TVib_wall = Twall
    refGrad()       = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}

void Foam::smoluchowskiJumpTVib::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);

    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("psi", "thermo2T:psi", psiName_);
    os.writeEntryIfDifferent<word>("mu", "thermo2T:mu", muName_);

    os.writeEntry("accommodationCoeff", accommodationCoeffVib_);
    TVibWall_.writeEntry("Twall", os);

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
        smoluchowskiJumpTVib
    );
}

// ************************************************************************* //