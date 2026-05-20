/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Class
    Foam::maxwellSlipU

Description
    Maxwell velocity slip boundary condition for the two-temperature
    non-equilibrium framework in nexaFoam.

    Applies the first-order Maxwell slip condition (Maxwell, 1879) to the
    velocity field U, acting through the partialSlip (Robin) base class.
    Two optional correction terms are included and independently switchable:

      - Thermal creep: tangential temperature gradient along the wall drives
        a slip velocity from cold toward hot regions, independent of shear.
        Uses the translational-rotational temperature TTR via TName_.

      - Curvature: non-equilibrium stress tensor correction for curved walls.
        Uses the Maxwell stress tensor tauMC, which must be computed and
        registered in the object registry by the solver before updateCoeffs()
        is called.

    The Robin blending formula from partialSlipFvPatchVectorField is:

        U_patch = f * Uwall + (1 - f) * U_interior  (tangential component)
        U_patch_normal = 0                           (always enforced)

    where:

        f   = 1 / (1 + delta * C1 * nu)
        C1  = sqrt(psi*pi/2) * (2 - alpha) / alpha
        nu  = mu / rho

    In the no-slip limit (Kn -> 0, C1 -> 0): f -> 1, U_patch -> Uwall.
    In the fully rarefied limit (Kn -> inf): f -> 0, U drifts to interior.

    The thermal creep correction is added to refValue (Uwall) as:

        refValue -= 3*nu/(4*TTR) * tangentialComponent(grad(TTR))

    The curvature correction is added to refValue as:

        refValue -= C1/rho * tangentialComponent(n . tauMC)

Usage
    Example boundary condition specification in 0/U:
    \verbatim
    wall
    {
        type                maxwellSlipU;
        accommodationCoeff  0.9;
        Uwall               uniform (0 0 0);
        thermalCreep        true;
        curvature           true;
        value               uniform (0 0 0);
    }
    \endverbatim

Note
    TName_ defaults to "TTR" — the translational-rotational temperature field
    as registered in nexaFoam. Override with T TTR; in the dict if different.
    tauMC must be computed by the solver before this BC fires — it is
    evaluated from muLam*dev2(T(fvc::grad(U))) in nexaFoam.C.

SourceFiles
    maxwellSlipU.C

\*---------------------------------------------------------------------------*/

#ifndef maxwellSlipU_H
#define maxwellSlipU_H

#include "partialSlipFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class maxwellSlipU Declaration
\*---------------------------------------------------------------------------*/

class maxwellSlipU
:
    public partialSlipFvPatchVectorField
{
    // Private data

        //- Translational-rotational temperature field name, default = "TTR"
        //  Used by the thermal creep term only. Must match the field name
        //  registered in the object registry by the thermo stack.
        word TName_;

        //- Density field name, default = "rho"
        word rhoName_;

        //- Compressibility field name, default = "thermo2T:psi"
        word psiName_;

        //- Dynamic viscosity field name, default = "thermo2T:mu"
        word muName_;

        //- Maxwell stress tensor field name, default = "tauMC"
        //  Must be registered by the solver before updateCoeffs() is called.
        word tauMCName_;

        //- Momentum accommodation coefficient (0, 1]
        scalar accommodationCoeff_;

        //- Prescribed wall velocity [m/s]
        //  Face-sized field to allow non-uniform moving walls.
        vectorField Uwall_;

        //- Include thermal creep term (default true)
        //  Switch off with thermalCreep false; in the BC dict.
        Switch thermalCreep_;

        //- Include surface curvature correction term (default true)
        //  Switch off with curvature false; in the BC dict.
        Switch curvature_;


public:

    //- Runtime type information
    TypeName("maxwellSlipU");


    // Constructors

        //- Construct from patch and internal field.
        //  Called by the framework for blank programmatic construction.
        //  Sets thermalCreep and curvature on, accommodationCoeff = 1.
        maxwellSlipU
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        //  Called by the run-time selection table when reading 0/U.
        //  All user parameters are read here: accommodationCoeff, Uwall,
        //  thermalCreep, curvature, and optional field name overrides.
        maxwellSlipU
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given maxwellSlipU onto a new patch.
        //  Called during AMR and dynamic mesh operations. Remaps all
        //  face-sized data (Uwall_) to the new patch face layout.
        maxwellSlipU
        (
            const maxwellSlipU&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference.
        //  Called by clone(iF) during decomposePar and parallel
        //  redistribution. All member variables including name strings
        //  must be copied explicitly.
        maxwellSlipU
        (
            const maxwellSlipU&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member functions

        // Mapping

            //- Map (and resize as needed) from self given a mapping object.
            //  Must remap Uwall_ to keep it consistent with new patch size.
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField.
            //  Must rmap Uwall_ to keep it consistent with source addressing.
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation

            //- Update the coefficients associated with the patch field.
            //  Computes C1, valueFraction, and refValue (Uwall + optional
            //  thermal creep and curvature corrections). All thermophysical
            //  quantities are derived from the thermo stack at runtime.
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //