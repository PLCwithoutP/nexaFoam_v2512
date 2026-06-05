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
    Foam::functionObjects::neStantonNumber

Group
    grpFieldFunctionObjects

Description
    Computes the Stanton number on wall patches from a wall heat-flux field
    (hy2Foam convention):

    \f[
        St = \frac{q_w}{\tfrac{1}{2}\,\rho_\infty\,|U_\infty|^3}
    \f]

    The denominator uses the freestream kinetic-energy flux as reference (it
    omits the total-enthalpy difference (H0 - Hw); for a cold hypersonic wall
    H0 ~ 1/2 U^2, so the two are close). The freestream reference is given via
    the dictionary; the wall heat-flux field is read from the registry at
    runtime or from disk under -postProcess.

Usage
    \verbatim
    StantonNumber
    {
        type            neStantonNumber;
        libs            ("libnexaFieldFunctionObjects.so");

        rhoInf          5.113e-4;
        magUInf         2764.5;

        wallHeatFlux    wallHeatFlux;   // must match neWallHeatFlux's field
        patches         (cone);
        result          StantonNumber;
        writeControl    writeTime;
        log             yes;
    }
    \endverbatim

SourceFiles
    neStantonNumber.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neStantonNumber_H
#define functionObjects_neStantonNumber_H

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
                       Class neStantonNumber Declaration
\*---------------------------------------------------------------------------*/

class neStantonNumber
:
    public fvMeshFunctionObject,
    public writeFile
{
    // Private Data

        //- Freestream (reference) density [kg/m^3]
        scalar rhoInf_;

        //- Freestream (reference) velocity magnitude [m/s]
        scalar magUInf_;

        //- Name of the wall heat-flux field to read
        word qName_;

        //- Name of the result field stored on the registry
        word resultName_;

        //- Wall patches to process (resolved in read(); empty => all walls)
        labelHashSet patchSet_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Return the wall heat-flux field from the registry, or read it from
        //  the current time directory if absent (post-process mode).
        tmp<volScalarField> lookupOrReadScalar(const word& fieldName) const;

        //- No copy construct / assignment
        neStantonNumber(const neStantonNumber&) = delete;
        void operator=(const neStantonNumber&) = delete;


public:

    //- Runtime type information
    TypeName("neStantonNumber");


    // Constructors

        //- Construct from Time and dictionary
        neStantonNumber
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~neStantonNumber() = default;


    // Member Functions

        //- Read the neStantonNumber data
        virtual bool read(const dictionary& dict);

        //- Calculate the Stanton-number field
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