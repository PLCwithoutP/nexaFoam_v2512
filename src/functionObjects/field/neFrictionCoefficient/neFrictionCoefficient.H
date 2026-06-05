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
    Foam::functionObjects::neFrictionCoefficient

Group
    grpFieldFunctionObjects

Description
    Computes the skin-friction coefficient on wall patches from a previously
    computed wall shear-stress (traction) field:

    \f[
        C_f = \frac{|\tau_{w,t}|}{\tfrac{1}{2}\,\rho_\infty\,|U_\infty|^2}
    \f]

    where \f$\tau_{w,t} = \tau_w - (\tau_w\cdot\hat{n})\hat{n}\f$ is the
    tangential part of the traction (the normal component is removed), and
    \f$\hat{n} = -S_f/|S_f|\f$.

    The freestream reference (rhoInf, magUInf) is user-specified via the
    dictionary, independent of any patch BC. The wall shear-stress field is
    read from the registry at runtime or from disk under -postProcess, so this
    object works in all three running modes.

Usage
    \verbatim
    frictionCoefficient
    {
        type            neFrictionCoefficient;
        libs            ("libnexaFieldFunctionObjects.so");

        rhoInf          5.113e-4;
        magUInf         2764.5;

        wallShearStress wallShearStress; // optional, default name shown
        patches         (cone);          // optional, default = all walls
        result          frictionCoefficient; // optional
        writeControl    writeTime;
        log             yes;
    }
    \endverbatim

SourceFiles
    neFrictionCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neFrictionCoefficient_H
#define functionObjects_neFrictionCoefficient_H

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
                   Class neFrictionCoefficient Declaration
\*---------------------------------------------------------------------------*/

class neFrictionCoefficient
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Freestream (reference) density [kg/m^3]
        scalar rhoInf_;

        //- Freestream (reference) velocity magnitude [m/s]
        scalar magUInf_;

        //- Name of the wall shear-stress field to read
        word wssName_;

        //- Name of the result field stored on the registry
        word resultName_;

        //- Wall patches to process (resolved in read(); empty => all walls)
        labelHashSet patchSet_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Return the wall shear-stress vector field from the registry, or
        //  read it from the current time directory if absent (post-process).
        tmp<volVectorField> lookupOrReadVector(const word& fieldName) const;

        //- No copy construct
        neFrictionCoefficient(const neFrictionCoefficient&) = delete;

        //- No copy assignment
        void operator=(const neFrictionCoefficient&) = delete;


public:

    //- Runtime type information
    TypeName("neFrictionCoefficient");


    // Constructors

        //- Construct from Time and dictionary
        neFrictionCoefficient
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~neFrictionCoefficient() = default;


    // Member Functions

        //- Read the neFrictionCoefficient data
        virtual bool read(const dictionary& dict);

        //- Calculate the friction-coefficient field
        virtual bool execute();

        //- Write the field and patch min/max/average
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //