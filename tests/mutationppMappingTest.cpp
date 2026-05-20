// =============================================================
//  Interface Mapping Test: N2+N | 2T | Thermal NonEq | Chem Eq
//
//  Maps every virtual function from:
//    - tcLibraryInterface.H
//    - compositionInterface.H
//  to its Mutation++ equivalent.
//
//  Legend in output:
//    [OK]  — confirmed Mutation++ function, directly mapped
//    [C]   — computed from confirmed functions (no direct call)
//    [TBD] — needs Mutation++ Doxygen API check before adapter impl
//    [N/A] — lifecycle/framework function, not applicable standalone
// =============================================================

#include "mutation++.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace std;

static const double R_UNIV     = 8.314462618;
static const double N_AV       = 6.02214076e23;
static const double M_N2       = 28.014e-3;
static const double M_N        = 14.007e-3;
static const double THETA_V_N2 = 3395.0;   // K  (RRHO, well-known)
static const double THETA_V_N  = 0.0;      // atomic

// ─────────────────────────────────────────────────────────────
//  Print helpers
// ─────────────────────────────────────────────────────────────

void printBanner(const string& title)
{
    const int w = 72;
    cout << "\n" << string(w,'#') << "\n";
    cout << "  " << title << "\n";
    cout << string(w,'#') << "\n";
}

void printSection(const string& s)
{
    const int w = 72;
    cout << "\n" << string(w,'-') << "\n";
    cout << "  " << s << "\n";
    cout << string(w,'-') << "\n";
    cout << "  " << setw(28) << left << "Interface fn"
                 << setw(34) << left << "Mutation++ mapping"
                 << "Value\n";
    cout << "  " << string(70,'.') << "\n";
}

void row(const string& fn, const string& mpp, double val,
         const string& unit, const string& tag)
{
    cout << "  " << setw(28) << left << fn
         << setw(34) << left << mpp
         << scientific << setprecision(4) << val
         << "  " << unit << "  " << tag << "\n";
}

void rowStr(const string& fn, const string& mpp,
            const string& val, const string& tag)
{
    cout << "  " << setw(28) << left << fn
         << setw(34) << left << mpp
         << val << "  " << tag << "\n";
}

void rowTBD(const string& fn, const string& note)
{
    cout << "  " << setw(28) << left << fn
         << note << "  [TBD]\n";
}

void rowNA(const string& fn, const string& note)
{
    cout << "  " << setw(28) << left << fn
         << note << "  [N/A]\n";
}

// ─────────────────────────────────────────────────────────────
//  Harmonic oscillator vibrational energy [J/kg]
// ─────────────────────────────────────────────────────────────

double esVib(double M, double thetaV, double T_v)
{
    if (thetaV <= 0.0 || T_v <= 0.0) return 0.0;
    return (R_UNIV / M) * thetaV / (std::exp(thetaV / T_v) - 1.0);
}

// ─────────────────────────────────────────────────────────────
//  Section 1: tcLibraryInterface
// ─────────────────────────────────────────────────────────────

void testTCLibrary(Mixture& mix)
{
    printSection("tcLibraryInterface");

    vector<double> Tvec(mix.nEnergyEqns());
    mix.getTemperatures(Tvec.data());
    double T_tr = Tvec[0];
    double T_v  = Tvec[1];

    // Thermodynamic state
    row("p()",    "mix.P()",              mix.P(),        "Pa",    "[OK]");
    row("TTR()",  "getTemperatures()[0]", T_tr,           "K",     "[OK]");
    row("TVib()", "getTemperatures()[1]", T_v,            "K",     "[OK]");

    // psi = rho/P   (ideal gas: p = rho * R_mix * T  =>  psi = 1/(R_mix*T) = rho/p)
    row("psi()",  "density()/P()",        mix.density()/mix.P(), "s2/m2", "[C]");

    // h() in neTCLib is hTR (TR enthalpy only, NOT total).
    // mix.mixtureHMass() = H_TR + e_V in 2T.
    // So hTR = mixtureHMass() - e_V_mix.
    int iN2 = mix.speciesIndex("N2");
    int iN  = mix.speciesIndex("N");
    double eV_mix = mix.Y()[iN2] * esVib(M_N2, THETA_V_N2, T_v)
                  + mix.Y()[iN]  * esVib(M_N,  THETA_V_N,  T_v);
    double hTR = mix.mixtureHMass() - eV_mix;
    row("h()  [= hTR]",   "mixtureHMass() - eV_mix", hTR,    "J/kg", "[C]");

    // eVib(): vibrational energy per unit mass of mixture
    row("eVib()", "SUM Y_i*esVib_i (HO)", eV_mix, "J/kg", "[C]");

    row("rho()", "density()",            mix.density(),  "kg/m3", "[OK]");

    // Transport
    row("mu()",  "viscosity()",          mix.viscosity(), "Pa.s", "[OK]");

    vector<double> lambda(mix.nEnergyEqns());
    mix.frozenThermalConductivityVector(lambda.data());
    row("kappaTR()",  "frozenThermalCond()[0]", lambda[0], "W/m-K", "[OK]");
    row("kappaVib()", "frozenThermalCond()[1]", lambda[1], "W/m-K", "[OK]");

    // fickDiffusionCoeff(): neTCLib returns a scalar mixture-level Fick coeff.
    // Mutation++ diffusion is per-species (full matrix or mixture-averaged Dm[i]).
    // Adapter must choose an averaging strategy.
    rowTBD("fickDiffusionCoeff()", "diffusionMatrix() diagonal — needs averaging");

    // gamma = Cp/Cv = Cp/(Cp - R_mix)   ideal gas relation
    double R_mix = mix.P() / (mix.density() * T_tr);
    double Cp    = mix.mixtureFrozenCpMass();
    row("gamma()", "CpMass/(CpMass-R_mix)", Cp/(Cp-R_mix), "-", "[C]");

    // Lifecycle — solver-framework operations, not standalone Mutation++ calls
    rowNA("correct()",            "solver-lifecycle (mesh/field update)");
    rowNA("correctVibEnergy()",   "solver-lifecycle (eVib field update)");
    rowNA("correctTVib()",        "solver-lifecycle (TVib field update)");
    rowNA("correctVibSource()",   "neTCLib Park VT  — TBD in adapter");
    rowNA("correctVibVibSource()","neTCLib Knab VV  — TBD in adapter");
    rowNA("correctCVSource()",    "neTCLib C-V coupling — TBD in adapter");

    // Configuration flags
    rowStr("twoTemperature()", "nEnergyEqns() > 1",
           mix.nEnergyEqns() > 1 ? "true" : "false", "[OK]");
    rowStr("chemistry()",      "nReactions() > 0",
           mix.nReactions()  > 0 ? "true" : "false", "[OK]");
    rowStr("chemistryCapable()","nReactions() > 0",
           mix.nReactions()  > 0 ? "true" : "false", "[OK]");

    // These three come from the OpenFOAM dict, not from Mutation++ runtime
    rowNA("speciesListFile()",   "from thermophysicalProperties dict");
    rowNA("reactionsListFile()", "from thermophysicalProperties dict");
    rowNA("chemistryReader()",   "from thermophysicalProperties dict");
}

// ─────────────────────────────────────────────────────────────
//  Section 2: compositionInterface
// ─────────────────────────────────────────────────────────────

void testComposition(Mixture& mix)
{
    printSection("compositionInterface");

    vector<double> Tvec(mix.nEnergyEqns());
    mix.getTemperatures(Tvec.data());
    double T_tr = Tvec[0];
    double T_v  = Tvec[1];
    double P    = mix.P();
    int    ns   = mix.nSpecies();
    int    iN2  = mix.speciesIndex("N2");
    int    iN   = mix.speciesIndex("N");

    // Per-species lookup arrays for this mixture
    vector<double> Mw(ns), thetaV(ns);
    Mw[iN2]     = M_N2;  thetaV[iN2] = THETA_V_N2;
    Mw[iN]      = M_N;   thetaV[iN]  = THETA_V_N;

    // ── Y() and species() ────────────────────────────────────

    cout << "\n  Y() → mix.Y()[i]    species() → mix.speciesName(i):\n";
    for (int i = 0; i < ns; i++)
        cout << "    Y_" << setw(5) << left << mix.speciesName(i)
             << ": " << scientific << setprecision(4) << mix.Y()[i]
             << "  [OK]\n";

    // ── Per-species two-column table ─────────────────────────

    auto pr2 = [&](const string& fn,
                   double vN2, double vN,
                   const string& note, const string& tag)
    {
        cout << "  " << setw(28) << left << fn
             << scientific << setprecision(3)
             << setw(14) << left << vN2
             << setw(14) << left << vN
             << note << "  " << tag << "\n";
    };

    cout << "\n";
    cout << "  " << setw(28) << left << "Function"
         << setw(14) << left << "N2"
         << setw(14) << left << "N"
         << "Mapping / Note\n";
    cout << "  " << string(70,'.') << "\n";

    // W(i) — molecular weight
    pr2("W(i) [kg/mol]",
        Mw[iN2], Mw[iN],
        "hardcoded — mix.M(i) or species obj needed", "[TBD]");

    // R(i) = R_UNIV / W(i)
    pr2("R(i) [J/kg-K]",
        R_UNIV/Mw[iN2], R_UNIV/Mw[iN],
        "R_UNIV / W(i)", "[C]");

    // thetaVib(i)
    pr2("thetaVib(i) [K]",
        thetaV[iN2], thetaV[iN],
        "RRHO DB — need species API", "[TBD]");

    // isSpecieMolecular(i)
    cout << "  " << setw(28) << left << "isSpecieMolecular(i)"
         << setw(14) << left << "true"
         << setw(14) << left << "false"
         << "from RRHO species type  [TBD]\n";

    // EsT(i, p, TTR) = (3/2)*(R/M)*T_tr   all species
    pr2("EsT(i) [J/kg]",
        (1.5)*(R_UNIV/Mw[iN2])*T_tr,
        (1.5)*(R_UNIV/Mw[iN]) *T_tr,
        "(3/2)*(R/M)*T_tr  (RRHO)", "[C]");

    // EsR(i, p, TTR): linear diatomic = R/M*T, monatomic = 0
    pr2("EsR(i) [J/kg]",
        (R_UNIV/Mw[iN2])*T_tr,  // N2: 2 rotational DOF
        0.0,                     // N:  monatomic
        "(R/M)*T_tr diatomic, 0 atomic  (RRHO)", "[C]");

    // EsVib(i, p, TTR, TVib, thetaVib) — harmonic oscillator
    pr2("EsVib(i) [J/kg]",
        esVib(Mw[iN2], thetaV[iN2], T_v),
        esVib(Mw[iN],  thetaV[iN],  T_v),
        "(R/M)*thetaV/(exp(thetaV/Tv)-1)", "[C]");

    // rho(i, p, TTR) = p*W_i/(R_UNIV*T_tr)   ideal gas
    pr2("rho(i) [kg/m3]",
        P*Mw[iN2]/(R_UNIV*T_tr),
        P*Mw[iN] /(R_UNIV*T_tr),
        "p*W_i/(R_UNIV*T_tr)  ideal gas", "[C]");

    // alphah(i) = kappaTR / (rho * Cp)
    // Mixture-level shown here — per-species needs species transport API
    {
        vector<double> lam(mix.nEnergyEqns());
        mix.frozenThermalConductivityVector(lam.data());
        double alpha = lam[0] / (mix.density() * mix.mixtureFrozenCpMass());
        cout << "  " << setw(28) << left << "alphah(i) [m2/s]"
             << "mix-level: " << scientific << setprecision(3) << alpha
             << "  (per-species needs species transport API)  [C/TBD]\n";
    }

    // Functions requiring Mutation++ species-level API (not in examples)
    cout << "\n  Remaining functions — need Mutation++ Doxygen API:\n";
    rowTBD("Hc(i)",       "formation enthalpy at 0K — in RRHO data");
    rowTBD("CpTR(i,p,T)", "per-species Cp — speciesCpOverR() or array fn");
    rowTBD("CvT(i,p,T)",  "CpTR - R/W(i)  once CpTR confirmed");
    rowTBD("H(i,p,T)",    "per-species H — speciesHOverRT() or array fn");
    rowTBD("Ha(i,p,T)",   "H(i) + Hc(i)");
    rowTBD("Hs(i,p,T)",   "H(i) - Hc(i)  sensible enthalpy");
    rowTBD("S(i,p,T)",    "per-species S — speciesSOverR() or array fn");
    rowTBD("G(i,p,T)",    "H - T*S  or speciesGOverRT()");
    rowTBD("mu(i,T)",     "species viscosity — Mutation++ transport module");
    rowTBD("kappaTR(i)",  "species TR conductivity — transport module");
    rowTBD("kappaVib(i)", "species V conductivity — transport module");
    rowTBD("TEVib(i,...)", "Newton inversion of HO: find Tv s.t. esVib(Tv)=given");
}

// ─────────────────────────────────────────────────────────────
//  main
// ─────────────────────────────────────────────────────────────

int main()
{
    printBanner(
        "Interface Mapping Test | N2+N | 2T | Thermal NonEq | Chemical Eq\n"
        "  tcLibraryInterface + compositionInterface → Mutation++ API");

    MixtureOptions opts;
    opts.setSpeciesDescriptor("N2 N");
    opts.setThermodynamicDatabase("RRHO");
    opts.setStateModel("ChemNonEqTTv");
    Mixture mix(opts);

    const double T_tr = 20000.0;
    const double T_v  = 300.0;
    const double n    = 1.0e22;    // [1/m3]
    const double xN2  = 0.8;
    const double xN   = 0.2;

    double rhoi[2];
    rhoi[mix.speciesIndex("N2")] = xN2 * n * M_N2 / N_AV;
    rhoi[mix.speciesIndex("N")]  = xN  * n * M_N  / N_AV;

    double T[2] = {T_tr, T_v};
    mix.setState(rhoi, T, 1);

    cout << "\n  State: T_tr=" << fixed << setprecision(0) << T_tr
         << " K  T_v=" << T_v << " K"
         << "  xN2=" << xN2 << "  xN=" << xN
         << "  n=" << scientific << n << " /m3\n";

    testTCLibrary(mix);
    testComposition(mix);

    const int w = 72;
    cout << "\n" << string(w,'=') << "\n";
    cout << "  KEY ARCHITECTURAL NOTE\n";
    cout << string(w,'=') << "\n";
    cout << "  neTCLib: each fn takes (speciei, p, TTR) — stateless,\n";
    cout << "           can be called at any (p,T) without side effects.\n\n";
    cout << "  Mutation++: STATEFUL — setState() must be called before\n";
    cout << "              every property query. One global Mixture object.\n\n";
    cout << "  => mutationPPAdapter.setState(cellI) must be called at the\n";
    cout << "     start of every cell loop before any property access.\n";
    cout << string(w,'=') << "\n";
    cout << "\n  [OK]=direct Mpp call  [C]=computed  [TBD]=needs Doxygen check"
         << "  [N/A]=framework\n\n";

    return 0;
}