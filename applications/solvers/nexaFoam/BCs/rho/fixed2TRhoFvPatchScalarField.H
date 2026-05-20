/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fixed2TRhoFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Foam::fixed2TRhoFvPatchScalarField

    This boundary condition provides a fixed density inlet condition for
    compressible solvers, where the density of calculated using:

    \f[
        \rho = \psi p
    \f]

    where
    \vartable
        p     | pressure [Pa]
        \rho  | density [kg/m3]
    \endvartable


Usage
    \table
        Property     | Description                | Required    | Default value
        p            | Pressure field name        | no          | p
        psi          | Compressibility field name | no          | thermo2T:psi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedRho;
    }
    \endverbatim

SourceFiles
    fixed2TRhoFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixed2TRhoFvPatchScalarField_H
#define fixed2TRhoFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixed2TRhoFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixed2TRhoFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

private:

    // Private data

        //- Pressure field name, default = "p"
        word pName_;

        //- Compressibility field name, default = "thermo2T:psi"
        word psiName_;


public:

    //- Runtime type information
    TypeName("fixed2TRho");


    // Constructors

        //- Construct from patch and internal field
        fixed2TRhoFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixed2TRhoFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixed2TRhoFvPatchScalarField
        // onto a new patch
        fixed2TRhoFvPatchScalarField
        (
            const fixed2TRhoFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixed2TRhoFvPatchScalarField
        (
            const fixed2TRhoFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        fixed2TRhoFvPatchScalarField
        (
            const fixed2TRhoFvPatchScalarField&,
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

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
