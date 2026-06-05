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
    Foam::functionObjects::neWallHeatFlux

Group
    grpFieldFunctionObjects

Description
    Computes the conductive wall heat flux [W/m^2] on wall patches:

    \f[
        q_w = \kappa_{TR}\,\frac{\partial T_{TR}}{\partial n}
            \;[ + \;\kappa_{Ve}\,\frac{\partial T_{Vib}}{\partial n} ]
    \f]

    Positive into the wall (cold-wall heating). The vibrational term is added
    when twoTemperature is true. Conductivities are read from written fields
    (publish them in the solver, like muLam); temperatures from TTR / TVib.
    Field lookups only, so the object works at runtime and under
    nexaFoam -postProcess / postProcess.

Usage
    \verbatim
    wallHeatFlux
    {
        type            neWallHeatFlux;
        libs            ("libnexaFieldFunctionObjects.so");

        kappaTR         kappaTR;    // optional, default "kappaTR"
        TTR             TTR;        // optional, default "TTR"

        // Optional: reference the wall-normal gradient to a fixed solid wall
        // temperature instead of the stored boundary face value. Use for
        // temperature-jump walls under -postProcess, whose written face value
        // collapses onto the near-wall cell (zero gas-side gradient on disk).
        // Twall           297.2;

        twoTemperature  false;      // optional, default false
        kappaVe         kappaVe;    // used only if twoTemperature
        TVib            TVib;       // used only if twoTemperature

        patches         (cone);     // optional, default = all wall patches
        result          wallHeatFlux; // optional, default = <function name>
        writeControl    writeTime;
        log             yes;
    }
    \endverbatim

SourceFiles
    neWallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neWallHeatFlux_H
#define functionObjects_neWallHeatFlux_H

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
                       Class neWallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class neWallHeatFlux
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Name of the translational-rotational conductivity field
        word kappaTRName_;

        //- Name of the translational-rotational temperature field
        word TTRName_;

        //- Optional fixed solid wall temperature [K] used as the gradient
        //  reference (for temperature-jump walls under -postProcess)
        scalar Twall_;

        //- True when Twall was supplied in the dictionary
        bool haveTwall_;

        //- Add the vibrational conduction term
        bool twoTemperature_;

        //- Name of the vibrational conductivity field
        word kappaVeName_;

        //- Name of the vibrational temperature field
        word TVibName_;

        //- Name of the result field stored on the registry
        word resultName_;

        //- Wall patches to process (resolved in read(); empty => all walls)
        labelHashSet patchSet_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Return a registered scalar field, or read it from the current time
        //  directory if absent (post-process mode).
        tmp<volScalarField> lookupOrReadScalar(const word& fieldName) const;

        //- No copy construct / assignment
        neWallHeatFlux(const neWallHeatFlux&) = delete;
        void operator=(const neWallHeatFlux&) = delete;


public:

    //- Runtime type information
    TypeName("neWallHeatFlux");


    // Constructors

        //- Construct from Time and dictionary
        neWallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~neWallHeatFlux() = default;


    // Member Functions

        //- Read the neWallHeatFlux data
        virtual bool read(const dictionary& dict);

        //- Calculate the wall heat-flux field
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