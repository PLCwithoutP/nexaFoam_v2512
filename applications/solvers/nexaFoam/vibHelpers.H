
/*---------------------------------------------------------------------------*\
    nexaFoam: vibrational-state helpers (nexaFoam Review Project, Phase 3)

    G1  exact inverse of the harmonic-oscillator e_v (replaces the floored
        library Newton, K1/A6)
    G2  energy-equivalent mixture vibrational temperature (replaces the
        mass-weighted mean, K4)
    G4  trace-molecule threshold

    The model is the library's: e_v = R theta/(exp(theta/T_v) - 1) with
    R = CpMCv of the perfect gas (HtoEthermo2T.H, EsVib). createFields.H
    checks at start-up that TvHO inverts the library EsVib for every molecule.
\*---------------------------------------------------------------------------*/

#ifndef vibHelpers_H
#define vibHelpers_H

#include <cmath>
#include "scalar.H"
#include "UList.H"

namespace Foam
{
namespace vibHelpers
{

//- A molecule with Y below this is a trace molecule (G4):
//  e_v = eps/(rho*Y) is not used, its T_v is set from the mixture T_v
constexpr scalar YTrace = 1e-10;


//- Harmonic-oscillator vibrational energy [J/kg]
inline scalar evHO(const scalar Tv, const scalar R, const scalar theta)
{
    if (Tv <= SMALL) return 0;
    const scalar x = theta/Tv;
    if (x > 700) return 0;
    return R*theta/std::expm1(x);
}


//- d(e_v)/d(T_v) of the harmonic oscillator [J/kg/K]
inline scalar cvvHO(const scalar Tv, const scalar R, const scalar theta)
{
    if (Tv <= SMALL) return 0;
    const scalar x = theta/Tv;
    if (x > 700) return 0;
    const scalar h = 0.5*x;
    const scalar r = h/std::sinh(h);
    return R*r*r;
}


//- G1: exact inverse T_v(e_v). Caller guarantees ev > VSMALL.
//  log1p keeps full accuracy in hot gas (R theta/ev -> 0).
inline scalar TvHO(const scalar ev, const scalar R, const scalar theta)
{
    return theta/std::log1p(R*theta/ev);
}


//- G2: energy-equivalent mixture T_v over the first n entries:
//  sum_i Y_i e_v,i(T) = sum_i Y_i ev_i.
//  Newton started at max_i Tv_i: f is increasing and convex and the root lies
//  in [min Tv_i, max Tv_i], so the iteration decreases monotonically to it.
//  maxRelErr is raised to the relative energy mismatch at the result.
inline scalar TvMix
(
    const label n,
    const UList<scalar>& Y,
    const UList<scalar>& R,
    const UList<scalar>& theta,
    const UList<scalar>& Tv,
    const UList<scalar>& ev,
    const scalar Tfallback,
    scalar& maxRelErr
)
{
    if (n == 0) return Tfallback;

    scalar Tmax = Tv[0];
    scalar Tmin = Tv[0];
    scalar target = 0;
    for (label i = 0; i < n; ++i)
    {
        Tmax = max(Tmax, Tv[i]);
        Tmin = min(Tmin, Tv[i]);
        target += Y[i]*ev[i];
    }

    // One molecule, or all at the same T_v: exact
    if (Tmax - Tmin <= 1e-12*Tmax) return Tmax;

    scalar T = Tmax;
    for (label iter = 0; iter < 100; ++iter)
    {
        scalar f = -target;
        scalar df = 0;
        for (label i = 0; i < n; ++i)
        {
            f  += Y[i]*evHO(T, R[i], theta[i]);
            df += Y[i]*cvvHO(T, R[i], theta[i]);
        }
        if (df <= VSMALL) break;        // vibrationally frozen point

        const scalar dT = f/df;
        T -= dT;
        if (mag(dT) <= 1e-12*T) break;
    }
    T = max(min(T, Tmax), Tmin);

    if (target > VSMALL)
    {
        scalar fEnd = -target;
        for (label i = 0; i < n; ++i)
        {
            fEnd += Y[i]*evHO(T, R[i], theta[i]);
        }
        maxRelErr = max(maxRelErr, mag(fEnd)/target);
    }

    return T;
}

} // End namespace vibHelpers
} // End namespace Foam

#endif