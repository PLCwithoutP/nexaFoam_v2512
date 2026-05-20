# nexaFoam
_non-equilibrium external aerodynamics_

---

**nexaFoam** is a density-based compressible flow solver specialized for **external hypersonic aerodynamics** under **thermal non-equilibrium**. It based upon `rhoCentralFoam` (OpenFOAM-v2512) and links against an in-house **neTCLib** (non-equilibrium ThermoChemical Library) for temperature-dependent thermodynamic/transport/kinetic data. 
Key capabilities:
- Two-temperature (**T** for translational–rotational, **Tᵥ** for vibrational) model 
- **No electrons, radicals or ions** (non-plasma, molecular nonequilibrium) 
- **Finite-rate reacting flows** (molecular chemistry) 

---

- **Physics**
  - Two-temperature formulation without electronic energy and without charged species
  - Vibrational–translational energy coupling via relaxation source terms
  - Temperature-dependent properties (cp(T), h(T), μ(T), κ(T), etc.) provided by **neTCLib**
- **Numerics**
  - Follows `rhoCentralFoam`’s density-based, central-upwind family of fluxes suitable for shocks
  - Spatial reconstruction with limiters (per standard OpenFOAM `fvSchemes`)
  - Time integration through OpenFOAM `ddtSchemes` (explicit/implicit choices supported by v2412)
- **Use cases**
  - Blunt/streamlined hypersonic bodies, wedges/ramps, fins
  - Shock layers with strong vibrational excitation
  - High-enthalpy entries where chemistry is important but **plasma effects are out of scope**


---

```
nexaFoam/
├─ applications/
│  └─ solvers/
│     └─ compressible/
│        └─ nexaFoam/           # solver sources (based on rhoCentralFoam)
├─ src/
│   └─ thermophysicalModels/
│      └─ neTCLib/              # in-house non-equilibrium thermo-chem library
│      └─ TCLibraryBase/        # thermo-chem library interface
├─ tutorials/                   # example/tutorial cases (optional)
├─ tests/                       # theoretical checks
├─ thirdParty/
│   └─ Mutationpp/
├─ README.md
```

---

## Requirements

- **OpenFOAM-v2512** (the solver relies on headers/ABIs of this release)
- **C++17** toolchain compatible with v2512
- Linux environment

---

## Limitations

- No electrons/ions, no plasma transport/conduction, no electron energy equation
- Reaction mechanisms restricted to **neutral** species
- Real-gas effects only if your neTCLib implementation supplies them
- Radiation and ablation not modeled unless explicitly added

---

## Acknowledgements

- Built on **OpenFOAM-v2512**’s `rhoCentralFoam`.
- neTCLib: in-house non-equilibrium thermochemical library for molecular gases.

