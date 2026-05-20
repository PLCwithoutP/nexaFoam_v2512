/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

Description
    FickModel: black-box helper to compute Fick diffusion mass fluxes.

    - Returns phiJ = J · Sf  [kg/s] for a given species mass fraction field Ys.
    - Also provides JsN = J · n [kg/(m^2 s)].

\*---------------------------------------------------------------------------*/

#ifndef FickModel_H
#define FickModel_H

#include "fvCFD.H"

namespace Foam
{

class FickModel
{
public:
    FickModel() = default;

    // J·n on faces (kg/(m^2 s))
    static inline tmp<surfaceScalarField> JsN
    (
        const volScalarField& rho,
        const volScalarField& Deff,
        const volScalarField& Ys
    )
    {
        // JsN = - (rho*Deff)_f * snGrad(Ys)
        return - fvc::interpolate(rho*Deff) * fvc::snGrad(Ys);
    }

    // J·Sf on faces (kg/s)
    static inline tmp<surfaceScalarField> phiJ
    (
        const volScalarField& rho,
        const volScalarField& Deff,
        const volScalarField& Ys
    )
    {
        const fvMesh& mesh = Ys.mesh();

        tmp<surfaceScalarField> tPhiJ
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phiJ_" + Ys.name(),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                - fvc::interpolate(rho*Deff) * fvc::snGrad(Ys) * mesh.magSf()
            )
        );

        // IMPORTANT: mark as oriented flux so it can be added to other oriented fluxes
        tPhiJ.ref().setOriented(true);

        return tPhiJ;
    }
};

} // namespace Foam

#endif

