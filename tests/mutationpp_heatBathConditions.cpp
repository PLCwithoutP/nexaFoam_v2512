// =============================================================
//  Mutation++ Heat Bath Capability Test
//  Covers Vatansever heat bath case matrix (5 cases)
//
//  Compile with mutationppTestRun.sh
// =============================================================

#include "mutation++.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace std;

// ─────────────────────────────────────────────────────────────
//  Print helpers
// ─────────────────────────────────────────────────────────────

void printCase(const string& title)
{
    const int w = 68;
    cout << "\n" << string(w, '=') << "\n";
    cout << "  " << title << "\n";
    cout << string(w, '=') << "\n";
}

void printSection(const string& title)
{
    cout << "\n  [" << title << "]\n";
}

void printModelInfo(Mixture& mix)
{
    printSection("Model");
    cout << "  nSpecies      : " << mix.nSpecies() << " -> ";
    for (int i = 0; i < mix.nSpecies(); i++)
        cout << mix.speciesName(i) << " ";
    cout << "\n";
    cout << "  nTemperatures : " << mix.nEnergyEqns() << "\n";
    cout << "  nReactions    : " << mix.nReactions() << "\n";
}

void printState(Mixture& mix)
{
    vector<double> Tvec(mix.nEnergyEqns());
    mix.getTemperatures(Tvec.data());

    printSection("State");
    cout << fixed << setprecision(2);
    cout << "  T_tr   [K]     : " << Tvec[0] << "\n";
    if (mix.nEnergyEqns() > 1)
        cout << "  T_v    [K]     : " << Tvec[1] << "\n";
    cout << "  P      [Pa]    : " << mix.P() << "\n";
    cout << scientific << setprecision(4);
    cout << "  rho    [kg/m3] : " << mix.density() << "\n";

    cout << "\n  Species mass fractions:\n";
    for (int i = 0; i < mix.nSpecies(); i++)
        cout << "    Y_" << setw(5) << left << mix.speciesName(i)
             << " : " << mix.Y()[i] << "\n";
}

void printThermo(Mixture& mix)
{
    printSection("Thermodynamics");
    cout << scientific << setprecision(4);
    cout << "  e_mix  [J/kg]   : " << mix.mixtureEnergyMass() << "\n";
    cout << "  h_mix  [J/kg]   : " << mix.mixtureHMass() << "\n";
    cout << "  Cp_fr  [J/kg-K] : " << mix.mixtureFrozenCpMass() << "\n";
}

void printTransport(Mixture& mix)
{
    printSection("Transport");
    cout << scientific << setprecision(4);
    cout << "  mu        [Pa.s]   : " << mix.viscosity() << "\n";

    vector<double> lambda(mix.nEnergyEqns());
    mix.frozenThermalConductivityVector(lambda.data());
    cout << "  lambda_TR [W/m-K]  : " << lambda[0] << "\n";
    if (mix.nEnergyEqns() > 1)
        cout << "  lambda_V  [W/m-K]  : " << lambda[1] << "\n";
}

void printChemistry(Mixture& mix)
{
    vector<double> wdot(mix.nSpecies(), 0.0);
    mix.netProductionRates(wdot.data());

    printSection("Chemistry source terms [kg/m3-s]");
    cout << scientific << setprecision(4);
    for (int i = 0; i < mix.nSpecies(); i++)
        cout << "    wdot_" << setw(5) << left << mix.speciesName(i)
             << " : " << wdot[i] << "\n";
}

// ─────────────────────────────────────────────────────────────
//  Shared constants
// ─────────────────────────────────────────────────────────────

const double N_A  = 6.02214076e23;  // Avogadro [1/mol]
const double M_N2 = 28.014e-3;      // kg/mol
const double M_N  = 14.007e-3;      // kg/mol
const double M_O2 = 31.999e-3;      // kg/mol

// ─────────────────────────────────────────────────────────────
//  Case 1: N2 only | 2T | Thermal NonEq | Chemical Eq (frozen)
// ─────────────────────────────────────────────────────────────
void case1_N2_2T()
{
    printCase("Case 1 | N2 | 2T | Thermal NonEq | Chemical Eq (frozen)");

    MixtureOptions opts;
    opts.setSpeciesDescriptor("N2");
    opts.setThermodynamicDatabase("RRHO");
    opts.setStateModel("ChemNonEqTTv");
    // No mechanism: composition is frozen. Pure 2T thermo + transport test.
    Mixture mix(opts);
    printModelInfo(mix);

    // Conditions: post-shock T_tr, frozen freestream T_v
    const double T_tr = 20000.0;   // K
    const double T_v  = 300.0;     // K
    const double n    = 1.0e22;    // number density [1/m3]

    double rhoi[1];
    rhoi[0] = n * M_N2 / N_A;      // pure N2

    double T[2] = {T_tr, T_v};
    mix.setState(rhoi, T, 1);

    printState(mix);
    printThermo(mix);
    printTransport(mix);
}

// ─────────────────────────────────────────────────────────────
//  Case 2: N2+N | 2T | Thermal NonEq | Chemical Eq (frozen)
// ─────────────────────────────────────────────────────────────
void case2_N2N_2T_thermalNonEq()
{
    printCase("Case 2 | N2+N | 2T | Thermal NonEq | Chemical Eq (frozen)");

    MixtureOptions opts;
    opts.setSpeciesDescriptor("N2 N");
    opts.setThermodynamicDatabase("RRHO");
    opts.setStateModel("ChemNonEqTTv");
    Mixture mix(opts);
    printModelInfo(mix);

    // 80% N2, 20% N by mole
    const double T_tr = 20000.0;
    const double T_v  = 300.0;
    const double n    = 1.0e22;
    const double xN2  = 0.8;
    const double xN   = 0.2;

    double rhoi[2];
    rhoi[mix.speciesIndex("N2")] = xN2 * n * M_N2 / N_A;
    rhoi[mix.speciesIndex("N")]  = xN  * n * M_N  / N_A;

    double T[2] = {T_tr, T_v};
    mix.setState(rhoi, T, 1);

    printState(mix);
    printThermo(mix);
    printTransport(mix);
}

// ─────────────────────────────────────────────────────────────
//  Case 3: N2+O2 | 2T | Thermal NonEq | Chemical Eq (frozen)
// ─────────────────────────────────────────────────────────────
void case3_N2O2_2T_thermalNonEq()
{
    printCase("Case 3 | N2+O2 | 2T | Thermal NonEq | Chemical Eq (frozen)");

    MixtureOptions opts;
    opts.setSpeciesDescriptor("N2 O2");
    opts.setThermodynamicDatabase("RRHO");
    opts.setStateModel("ChemNonEqTTv");
    Mixture mix(opts);
    printModelInfo(mix);

    // Air-like: 79% N2, 21% O2 by mole
    const double T_tr = 20000.0;
    const double T_v  = 300.0;
    const double n    = 1.0e22;
    const double xN2  = 0.79;
    const double xO2  = 0.21;

    double rhoi[2];
    rhoi[mix.speciesIndex("N2")] = xN2 * n * M_N2 / N_A;
    rhoi[mix.speciesIndex("O2")] = xO2 * n * M_O2 / N_A;

    double T[2] = {T_tr, T_v};
    mix.setState(rhoi, T, 1);

    printState(mix);
    printThermo(mix);
    printTransport(mix);
}

// ─────────────────────────────────────────────────────────────
//  Case 4: Air-5 | 1T | Chemical NonEq
//  Uses built-in air_5 mixture + mechanism
// ─────────────────────────────────────────────────────────────
void case4_Air5_1T_chemNonEq()
{
    printCase("Case 4 | Air-5 | 1T | Chemical NonEq");

    // air_5 already carries its mechanism file (5 reactions)
    MixtureOptions opts("air_5");
    opts.setThermodynamicDatabase("RRHO");
    opts.setStateModel("ChemNonEq1T");
    Mixture mix(opts);
    printModelInfo(mix);

    // Non-equilibrium initial state: N2+O2 only, no dissociation products yet
    // At T=5000K this will produce non-zero omega_dot showing N2/O2 dissociation
    const double T   = 5000.0;    // K
    const double n   = 1.0e22;    // number density [1/m3]
    const double xN2 = 0.79;
    const double xO2 = 0.21;

    // air_5 has 5 species: N O NO N2 O2
    // All set to zero first, then fill N2 and O2 only
    double rhoi[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    rhoi[mix.speciesIndex("N2")] = xN2 * n * M_N2 / N_A;
    rhoi[mix.speciesIndex("O2")] = xO2 * n * M_O2 / N_A;

    mix.setState(rhoi, &T, 1);

    printState(mix);
    printThermo(mix);
    printTransport(mix);
    printChemistry(mix);
}

// ─────────────────────────────────────────────────────────────
//  Case 5: N2+N | 2T | Thermochemical NonEq
// ─────────────────────────────────────────────────────────────
void case5_N2N_2T_thermoChemNonEq()
{
    printCase("Case 5 | N2+N | 2T | Thermochemical NonEq");

    MixtureOptions opts;
    opts.setSpeciesDescriptor("N2 N");
    opts.setThermodynamicDatabase("RRHO");
    opts.setStateModel("ChemNonEqTTv");
    // Mechanism: uncomment once you verify the name in:
    //   ls $MPP_DATA_DIRECTORY/mechanisms/
    // opts.setMechanism("N2_diss");
    Mixture mix(opts);
    printModelInfo(mix);

    // 90% N2, 10% N by mole, T_tr >> T_v
    const double T_tr = 20000.0;
    const double T_v  = 5000.0;
    const double n    = 1.0e22;
    const double xN2  = 0.9;
    const double xN   = 0.1;

    double rhoi[2];
    rhoi[mix.speciesIndex("N2")] = xN2 * n * M_N2 / N_A;
    rhoi[mix.speciesIndex("N")]  = xN  * n * M_N  / N_A;

    double T[2] = {T_tr, T_v};
    mix.setState(rhoi, T, 1);

    printState(mix);
    printThermo(mix);
    printTransport(mix);

    // Chemistry: uncomment after mechanism is set above
    // printChemistry(mix);

    cout << "\n  NOTE: mechanism not loaded — omega_dot skipped.\n";
    cout << "  Run: ls $MPP_DATA_DIRECTORY/mechanisms/ to find N2 mechanism.\n";
}

// ─────────────────────────────────────────────────────────────
//  main
// ─────────────────────────────────────────────────────────────
int main()
{
    const int w = 68;
    cout << "\n" << string(w, '#') << "\n";
    cout << "  Mutation++ Heat Bath Capability Test\n";
    cout << "  Vatansever heat bath case matrix\n";
    cout << string(w, '#') << "\n";

    case1_N2_2T();
    case2_N2N_2T_thermalNonEq();
    case3_N2O2_2T_thermalNonEq();
    case4_Air5_1T_chemNonEq();
    case5_N2N_2T_thermoChemNonEq();

    cout << "\n" << string(w, '#') << "\n";
    cout << "  All cases completed.\n";
    cout << string(w, '#') << "\n\n";

    return 0;
}