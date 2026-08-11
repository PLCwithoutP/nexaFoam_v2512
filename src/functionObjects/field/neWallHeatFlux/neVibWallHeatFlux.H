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
    Foam::functionObjects::neVibWallHeatFlux

Group
    grpFieldFunctionObjects

Description
    Computes the vibrational-electronic conductive wall heat flux [W/m^2]
    on wall patches, for two-temperature runs:

    \f[
        q_{w,\text{vib}} = \kappa_{Ve}\,\frac{\partial T_{Vib}}{\partial n}
    \f]

    Positive into the wall (cold-wall heating). Conductivity is read from a
    written field (publish it in the solver, like muLam); temperature from
    TVib. Field lookups only, so the object works at runtime and under
    nexaFoam -postProcess / postProcess.

    Unlike the translational-rotational term in neWallHeatFlux, kappaVe and
    TVib are not optional here: this object exists specifically to report
    the vibrational term, so a missing field is a configuration error and
    lookupOrReadScalar's MUST_READ will fail fatally rather than silently
    skipping the term. Only add this function object to 2T cases; for 1T /
    pure-gas / non-reacting runs, use neWallHeatFlux alone.

Usage
    \verbatim
    wallHeatFluxVib
    {
        type            neVibWallHeatFlux;
        libs            ("libnexaFieldFunctionObjects.so");

        kappaVe         kappaVe;    // vib conductivity field (default "kappaVe")
        TVib            TVib;       // vib temperature field (default "TVib")

        // Optional: reference the wall-normal gradient to a fixed solid wall
        // vibrational temperature instead of the stored boundary face value.
        // Use for temperature-jump walls under -postProcess, whose written
        // face value collapses onto the near-wall cell (zero gas-side
        // gradient on disk).
        // TVibWall        297.2;

        patches         (cone);        // optional, default = all wall patches
        result          wallHeatFluxVib; // optional, default = <function name>
        writeControl    writeTime;
        log             yes;
    }
    \endverbatim

SourceFiles
    neVibWallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neVibWallHeatFlux_H
#define functionObjects_neVibWallHeatFlux_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFields.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class neVibWallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class neVibWallHeatFlux
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Name of the vibrational-electronic conductivity field
        word kappaVeName_;

        //- Name of the vibrational temperature field
        word TVibName_;

        //- Optional fixed solid wall vibrational temperature [K] used as the
        //  gradient reference (for temperature-jump walls under -postProcess)
        scalar TVibWall_;

        //- True when TVibWall was supplied in the dictionary
        bool haveTVibWall_;

        //- Name of the result field stored on the registry
        word resultName_;

        //- Wall patches to process (resolved in read(); empty => all walls)
        labelHashSet patchSet_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Return a registered scalar field, or read it from the current time
        //  directory if absent (post-process mode). Fatal if not found in
        //  either place, since kappaVe/TVib are mandatory for this object.
        tmp<volScalarField> lookupOrReadScalar(const word& fieldName) const;

        //- Accumulate kappa*(T_cell - T_face)/delta into the wall boundary of
        //  the result field. When useWallT is true, wallT is used as the
        //  surface reference instead of the stored boundary face value.
        void addConductionTerm
        (
            volScalarField::Boundary& qBf,
            const volScalarField& kappa,
            const volScalarField& T,
            const bool useWallT,
            const scalar wallT
        ) const;

        //- No copy construct / assignment
        neVibWallHeatFlux(const neVibWallHeatFlux&) = delete;
        void operator=(const neVibWallHeatFlux&) = delete;


public:

    //- Runtime type information
    TypeName("neVibWallHeatFlux");


    // Constructors

        //- Construct from Time and dictionary
        neVibWallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~neVibWallHeatFlux() = default;


    // Member Functions

        //- Read the neVibWallHeatFlux data
        virtual bool read(const dictionary& dict);

        //- Calculate the vibrational wall heat-flux field
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
