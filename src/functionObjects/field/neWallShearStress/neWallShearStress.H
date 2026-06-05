/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    nexaFoam field function object.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::functionObjects::neWallShearStress

Group
    grpFieldFunctionObjects

Description
    Computes the laminar wall shear-stress (traction) vector on wall patches:

    \f[
        \tau_w = \hat{n}_\mathrm{in} \cdot R_\mathrm{eff}, \qquad
        R_\mathrm{eff} = -\mu\,\mathrm{dev}(2\,\mathrm{symm}(\nabla U))
    \f]

    where \f$\hat{n}_\mathrm{in} = -S_f/|S_f|\f$ is the inward wall normal.
    \f$R_\mathrm{eff}\f$ is the laminar deviatoric viscous stress (equal to the
    standard laminar devRhoReff); turbulence is not used.

    The dynamic viscosity is read from the registered/written field 'muLam'
    (publish it in the solver), and U from the registered/written field 'U'.
    Only field lookups are used, so the object works at runtime and under
    nexaFoam -postProcess / postProcess.

Usage
    \verbatim
    wallShearStress
    {
        type            neWallShearStress;
        libs            ("libnexaFieldFunctionObjects.so");

        mu              muLam;          // optional, default "muLam"
        U               U;             // optional, default "U"
        patches         (cone);        // optional, default = all wall patches
        result          wallShearStress; // optional, default = <function name>
        writeControl    writeTime;
        log             yes;
    }
    \endverbatim

SourceFiles
    neWallShearStress.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neWallShearStress_H
#define functionObjects_neWallShearStress_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class neWallShearStress Declaration
\*---------------------------------------------------------------------------*/

class neWallShearStress
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Name of the dynamic-viscosity field
        word muName_;

        //- Name of the velocity field
        word UName_;

        //- Name of the result field stored on the registry
        word resultName_;

        //- Wall patches to process (resolved in read(); empty => all walls)
        labelHashSet patchSet_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Return a registered scalar field, or read it from the current time
        //  directory if absent (post-process mode). Non-owning tmp at runtime,
        //  owning tmp offline.
        tmp<volScalarField> lookupOrReadScalar(const word& fieldName) const;

        //- As above, for a vector field (used for U).
        tmp<volVectorField> lookupOrReadVector(const word& fieldName) const;

        //- No copy construct
        neWallShearStress(const neWallShearStress&) = delete;

        //- No copy assignment
        void operator=(const neWallShearStress&) = delete;


public:

    //- Runtime type information
    TypeName("neWallShearStress");


    // Constructors

        //- Construct from Time and dictionary
        neWallShearStress
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~neWallShearStress() = default;


    // Member Functions

        //- Read the neWallShearStress data
        virtual bool read(const dictionary& dict);

        //- Calculate the wall shear-stress field
        virtual bool execute();

        //- Write the field and patch magnitude min/max/average
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //