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
    Computes the translational-rotational conductive wall heat flux [W/m^2]
    on wall patches:

    \f[
        q_w = \kappa_{TR}\,\frac{\partial T_{TR}}{\partial n}
    \f]

    Positive into the wall (cold-wall heating). Conductivity is read from a
    written field (publish it in the solver, like muLam); temperature from
    TTR. Field lookups only, so the object works at runtime and under
    nexaFoam -postProcess / postProcess.

    For the vibrational-electronic counterpart (kappaVe / TVib), see
    neVibWallHeatFlux.

Usage
    \verbatim
    wallHeatFlux
    {
        type            neWallHeatFlux;
        libs            ("libnexaFieldFunctionObjects.so");
        qFields         (qWallTR qWallVib);    // 1T: (qWallTR)
        writeControl    writeTime;
    }
    \endverbatim

SourceFiles
    neWallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neWallHeatFlux_H
#define functionObjects_neWallHeatFlux_H

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

        //- Name of the result field stored on the registry
        word resultName_;

        //- Wall patches to process (resolved in read(); empty => all walls)
        labelHashSet patchSet_;

        //- Solver wall heat-flux fields to report (P4-08); empty: kappa dT/dn
        wordList qFieldNames_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Return a registered scalar field, or read it from the current time
        //  directory if absent (post-process mode).
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
