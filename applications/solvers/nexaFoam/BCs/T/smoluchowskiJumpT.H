/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::smoluchowskiJumpT

Description
    Smoluchowski temperature jump boundary condition for the two-temperature
    non-equilibrium framework in nexaFoam.

    Applies the temperature jump condition (Smoluchowski, 1898) in enthalpy
    space, acting on the hTR field rather than TTR directly. TTR is a
    diagnostic field owned by the thermo stack and recomputed by
    thermo2T.correct() — imposing a BC on TTR directly would be overwritten
    on every thermo update. The jump is therefore encoded into hTR at the
    boundary so that thermo inversion propagates the physically correct
    jumped TTR into the interior.

    The Robin (mixed) blending formula is:

        hTR_patch = f * hTRwall + (1 - f) * (hTR_internal + refGrad/delta)

    where:

        f       = 1 / (1 + delta * C2)
        C2      = (mu/rho) * sqrt(psi*pi/2) * 2*gamma/((gamma+1)*Pr)
                  * (2 - alpha) / alpha
        hTRwall = Cp * Twall

    Pr and gamma are computed face-by-face from the thermo stack — no
    fixed scalar assumptions are made. This is required for thermally
    perfect non-equilibrium multispecies mixtures where both quantities
    vary with temperature and composition.

    In the no-jump limit (Kn -> 0, C2 -> 0): f -> 1, hTR_patch -> hTRwall.
    In the fully rarefied limit (Kn -> inf): f -> 0, hTR drifts to interior.

Usage
    Example boundary condition specification in 0/hTR:
    \verbatim
    wall
    {
        type                smoluchowskiJumpT;
        accommodationCoeff  0.9;
        Twall               uniform 300.0;
        value               uniform 0.0;
    }
    \endverbatim

Note
    gamma and Pr are not user inputs — they are derived from the thermo
    stack at every updateCoeffs() call. Do not specify them in the dict.

SourceFiles
    smoluchowskiJumpT.C

\*---------------------------------------------------------------------------*/

#ifndef smoluchowskiJumpT_H
#define smoluchowskiJumpT_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class smoluchowskiJumpT Declaration
\*---------------------------------------------------------------------------*/

class smoluchowskiJumpT
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Velocity field name, default = "U"
        word UName_;

        //- Density field name, default = "rho"
        word rhoName_;

        //- Compressibility field name, default = "thermo2T:psi"
        word psiName_;

        //- Dynamic viscosity field name, default = "thermo2T:mu"
        word muName_;

        //- Thermal accommodation coefficient (0, 1]
        scalar accommodationCoeff_;

        //- Prescribed wall surface temperature [K]
        //  Face-sized field to allow non-uniform wall temperature.
        scalarField Twall_;


public:

    //- Runtime type information
    TypeName("smoluchowskiJumpT");


    // Constructors

        //- Construct from patch and internal field.
        //  Called by the framework for blank programmatic construction.
        //  Sets safe zero defaults for refValue, refGrad, valueFraction.
        smoluchowskiJumpT
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary.
        //  Called by the run-time selection table when reading 0/hTR.
        //  All user parameters (accommodationCoeff, Twall) are read here.
        smoluchowskiJumpT
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given smoluchowskiJumpT onto a new patch.
        //  Called during AMR and dynamic mesh operations. Remaps all
        //  face-sized data (Twall_) to the new patch face layout.
        smoluchowskiJumpT
        (
            const smoluchowskiJumpT&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference.
        //  Called by clone(iF) during decomposePar and parallel
        //  redistribution. All member variables including name strings
        //  must be copied explicitly.
        smoluchowskiJumpT
        (
            const smoluchowskiJumpT&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        // Mapping

            //- Map (and resize as needed) from self given a mapping object.
            //  Must remap Twall_ to keep it consistent with new patch size.
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField.
            //  Must rmap Twall_ to keep it consistent with source addressing.
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );

            virtual bool fixesValue() const 
            { 
                return true; 
            }
        // Evaluation

            //- Update the coefficients associated with the patch field.
            //  Computes C2, valueFraction, and refValue (hTRwall = Cp*Twall)
            //  from the thermo stack. Pr and gamma are derived locally —
            //  no fixed scalar assumptions.
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //