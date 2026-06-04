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
    Foam::functionObjects::nePressureCoefficient

Group
    grpFieldFunctionObjects

Description
    Computes the surface pressure coefficient on wall patches:

    \f[
        C_p = \frac{p_\mathrm{wall}}{\tfrac{1}{2}\,\rho_\infty\,|U_\infty|^2}
    \f]

    The freestream reference is user-specified via the dictionary entries
    rhoInf and magUInf (it is NOT read from any patch, so it is independent
    of boundary-condition types). Only the registered field p is read, so the
    object works at runtime and under nexaFoam -postProcess / postProcess.

Usage
    \verbatim
    cp_surface
    {
        type            nePressureCoefficient;
        libs            ("libnexaFieldFunctionObjects.so");

        rhoInf          5.113e-4;   // freestream density [kg/m^3]   (required)
        magUInf         2764.5;     // freestream speed   [m/s]      (required)

        patches         (cone);     // optional; default = all wall patches
        result          cp_surface; // optional; default = <function name>
        writeControl    writeTime;
        log             yes;
    }
    \endverbatim

SourceFiles
    nePressureCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_nePressureCoefficient_H
#define functionObjects_nePressureCoefficient_H

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
                   Class nePressureCoefficient Declaration
\*---------------------------------------------------------------------------*/

class nePressureCoefficient
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Freestream (reference) density [kg/m^3]
        scalar rhoInf_;

        //- Freestream (reference) velocity magnitude [m/s]
        scalar magUInf_;

        //- Name of the result field stored on the registry
        word resultName_;

        //- Wall patches to process (resolved in read(); empty => all walls)
        labelHashSet patchSet_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Return the registered scalar field p, or read it from the current
        //  time directory if it is not in the registry (post-process mode).
        //  The returned tmp references the registry field at runtime and owns
        //  the read field offline; valid for the caller's scope either way.
        tmp<volScalarField> lookupOrReadScalar(const word& fieldName) const;

        //- No copy construct
        nePressureCoefficient(const nePressureCoefficient&) = delete;

        //- No copy assignment
        void operator=(const nePressureCoefficient&) = delete;


public:

    //- Runtime type information
    TypeName("nePressureCoefficient");


    // Constructors

        //- Construct from Time and dictionary
        nePressureCoefficient
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~nePressureCoefficient() = default;


    // Member Functions

        //- Read the nePressureCoefficient data
        virtual bool read(const dictionary& dict);

        //- Calculate the nePressureCoefficient field
        virtual bool execute();

        //- Write the nePressureCoefficient field and patch min/max/average
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //