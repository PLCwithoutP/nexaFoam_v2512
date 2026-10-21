/*---------------------------------------------------------------------------*\
  nexaFoam (OpenFOAM-v2512): ported from IHIS-v2312 (nexaFoam Review Project,
  Phase 4, D4.4/Q4.7). Only the access to the thermo state differs (marked P4).
  Derived from rhoCentralFoam/BCs/U/maxwellSlipUFvPatchVectorField
  (Copyright (C) 2011-2016 OpenFOAM Foundation, 2020 OpenCFD Ltd., GPL v3)
\*---------------------------------------------------------------------------*/

#include "neMaxwellSlipUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "neRarefiedWall.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neMaxwellSlipUFvPatchVectorField::neMaxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    TName_("TTR"),
    tauMCName_("tauMC"),
    accommodationCoeff_(1.0),
    omega_(0.74),       // only without a dictionary; the dictionary requires omega
    Uwall_(p.size(), Zero),
    thermalCreep_(false),
    curvature_(false)
{}


Foam::neMaxwellSlipUFvPatchVectorField::neMaxwellSlipUFvPatchVectorField
(
    const neMaxwellSlipUFvPatchVectorField& mspvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(mspvf, p, iF, mapper),
    TName_(mspvf.TName_),
    tauMCName_(mspvf.tauMCName_),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    omega_(mspvf.omega_),
    Uwall_(mspvf.Uwall_, mapper),       // mapped (stock copies it unmapped)
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


Foam::neMaxwellSlipUFvPatchVectorField::neMaxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    TName_(dict.getOrDefault<word>("T", "TTR")),
    tauMCName_(dict.getOrDefault<word>("tauMC", "tauMC")),
    accommodationCoeff_
    (
        neRarefiedWall::readAccommodation(dict, "accommodationCoeff")
    ),
    omega_(neRarefiedWall::readOmega(dict)),
    Uwall_(p.size(), Zero),
    thermalCreep_(dict.getOrDefault<Switch>("thermalCreep", false)),
    curvature_(dict.getOrDefault<Switch>("curvature", false))
{
    if (dict.found("Uwall"))
    {
        Uwall_ = vectorField("Uwall", dict, p.size());
    }

    // Until the coefficients are first available: U_b = U_w (no-slip)
    refValue() = Uwall_;
    valueFraction() = scalar(1);

    if (this->readValueEntry(dict))
    {
        // Restart: keep the written state
        const auto* hasRefValue = dict.findEntry("refValue", keyType::LITERAL);
        const auto* hasFrac = dict.findEntry("valueFraction", keyType::LITERAL);

        if (hasRefValue && hasFrac)
        {
            this->refValue().assign(*hasRefValue, p.size());
            this->valueFraction().assign(*hasFrac, p.size());
        }
    }
    else
    {
        fvPatchVectorField::operator=(Uwall_);
    }
}


Foam::neMaxwellSlipUFvPatchVectorField::neMaxwellSlipUFvPatchVectorField
(
    const neMaxwellSlipUFvPatchVectorField& mspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(mspvf, iF),
    TName_(mspvf.TName_),
    tauMCName_(mspvf.tauMCName_),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    omega_(mspvf.omega_),
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::neMaxwellSlipUFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    partialSlipFvPatchVectorField::autoMap(m);
    Uwall_.autoMap(m);
}


void Foam::neMaxwellSlipUFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    partialSlipFvPatchVectorField::rmap(ptf, addr);

    const auto& mptf = refCast<const neMaxwellSlipUFvPatchVectorField>(ptf);
    Uwall_.rmap(mptf.Uwall_, addr);
}


void Foam::neMaxwellSlipUFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!neRarefiedWall::ready(db()))
    {
        partialSlipFvPatchVectorField::updateCoeffs();
        return;
    }

    // P4: thermo state from the nexaFoam registry (p, psi, muLam patch values)
    const tcLibraryInterface& thermo = neRarefiedWall::thermo(db());

    const label patchi = patch().index();

    const scalarField& pp = patch().lookupPatchField<volScalarField>("p");
    const scalarField& ppsi = thermo.psi().boundaryField()[patchi];
    const scalarField& pmu = patch().lookupPatchField<volScalarField>("muLam");

    // C1 = lambda (2 - s)/s  [m]
    const scalarField C1
    (
        neRarefiedWall::meanFreePathVHS(pmu, pp, ppsi, omega_)
       *neRarefiedWall::accommodationFactor(accommodationCoeff_)
    );

    valueFraction() = 1.0/(1.0 + patch().deltaCoeffs()*C1);

    refValue() = Uwall_;

    if (thermalCreep_ || curvature_)
    {
        const vectorField n(patch().nf());

        if (thermalCreep_)
        {
            const volScalarField& vsfT =
                db().lookupObject<volScalarField>(TName_);
            const scalarField& pT = vsfT.boundaryField()[patchi];
            const vectorField gradpT(fvc::grad(vsfT)().boundaryField()[patchi]);

            // nu = mu/rho_w, rho_w = p psi
            const scalarField pnu(pmu/(pp*ppsi));

            // Physical sign: the gas creeps toward higher temperature
            // (Maxwell 1879; Kennard 1938). Stock OpenFOAM and hy2Foam use -=.
            refValue() += 0.75*pnu/pT*transform(I - n*n, gradpT);
        }

        if (curvature_)
        {
            const fvPatchTensorField& ptauMC =
                patch().lookupPatchField<volTensorField>(tauMCName_);

            // lambda (2 - s)/s / mu * tangential part of n.tauMC
            refValue() -= C1/pmu*transform(I - n*n, (n & ptauMC));
        }
    }

    partialSlipFvPatchVectorField::updateCoeffs();
}


void Foam::neMaxwellSlipUFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("T", "TTR", TName_);
    os.writeEntryIfDifferent<word>("tauMC", "tauMC", tauMCName_);

    os.writeEntry("accommodationCoeff", accommodationCoeff_);
    os.writeEntry("omega", omega_);
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
        neMaxwellSlipUFvPatchVectorField
    );
}

// ************************************************************************* //