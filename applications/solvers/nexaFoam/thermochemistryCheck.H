/*---------------------------------------------------------------------------*\
    thermochemistryCheck.H

    Reads thermophysicalProperties to establish two solver-wide booleans:

      mixtureCheck         — true if mixture is multi-component (per-species
                             property access required in composition interface)
      useVibVibRelaxation  — true if Knab VV relaxation source terms are active

    Selection logic branches on the 'library' keyword:

      library neTCLib (or default):
        mixtureCheck       ← thermoType.mixture keyword
                             ("multiComponent2TMixture" or "neReactingMixture"
                              → true; "pure2TMixture" → false)
        useVibVibRelaxation ← VVrelaxation keyword (default false)

      library mutationpp:
        mixtureCheck       ← explicit 'mixtureCheck' boolean in
                             thermophysicalProperties (default true).
                             Cannot be derived from nSpecies() here because
                             the adapter has not been constructed yet.
        useVibVibRelaxation ← always false (Knab VV not in Mutation++ API)
\*---------------------------------------------------------------------------*/

IOdictionary thermoDict_
(
    IOobject
    (
        "thermophysicalProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    )
);

// ── Read library keyword ──────────────────────────────────────────────────────

const word libName_ =
    thermoDict_.lookupOrDefault<word>("library", "neTCLib");
const bool inviscid = 
    thermoDict_.lookupOrDefault<bool>("inviscid", false);

bool mixtureCheck        = false;
bool useVibVibRelaxation = false;


// ── Branch on backend ─────────────────────────────────────────────────────────

if (libName_ == "mutationpp")
{
    // On the Mutation++ path, mixtureCheck controls which vibrational
    // energy equation branch is taken (vibEnergyEquation.H).
    // The !mixtureCheck branch uses QVzero — no VT source ever.
    // The mixtureCheck branch applies QVT per species.
    // Therefore: twoTemperature → mixtureCheck must always be true.
    // Single-species thermal nonequilibrium still uses the per-species
    // path (Y.N2 = 1 everywhere, one equation, VT source applied).

    const bool use2T =
        thermoDict_.lookupOrDefault<bool>("twoTemperature", false);

    const bool chemistryRequested =
        thermoDict_.lookupOrDefault<bool>("chemistry", false);
    
    mixtureCheck =   use2T
                || chemistryRequested
                || thermoDict_.lookupOrDefault<bool>("mixtureCheck", false);

    useVibVibRelaxation = false;

    Info<< "Mutation++ backend detected." << nl
        << "  mixtureCheck        = "
        << (mixtureCheck ? "true" : "false")
        << (use2T        ? " (forced true: twoTemperature active)"  :
            chemistryRequested ? " (forced true: chemistry active)"      : "") << nl
        << "  useVibVibRelaxation = false (VV not in Mutation++ API)" << nl
        << endl;

    Info<< "  energyMode = vibrational-electronic (Ve)"
        << " — Mutation++ RRHO includes electronic states" << nl;
}
else
{
    // ── neTCLib path (default) ────────────────────────────────────────────────
    //
    // mixtureCheck is derived from thermoType.mixture keyword.
    // useVibVibRelaxation is read from VVrelaxation keyword.

    if (thermoDict_.found("thermoType"))
    {
        const dictionary& thermoTypeDict = thermoDict_.subDict("thermoType");

        if (thermoTypeDict.found("mixture"))
        {
            const word mixtureType(thermoTypeDict.lookup("mixture"));

            if (mixtureType == "multiComponent2TMixture")
            {
                mixtureCheck = true;
                Info<< "Multiple species are found in mixture." << nl
                    << endl;
            }
            else if (mixtureType == "neReactingMixture")
            {
                mixtureCheck = true;
                Info<< "Reacting multi-specie mixture is detected." << nl
                    << endl;
            }
            else if (mixtureType == "pure2TMixture")
            {
                mixtureCheck = false;
                Info<< "Pure mixture is detected." << nl
                    << endl;
            }
            else
            {
                WarningInFunction
                    << "Unknown thermoType.mixture = " << mixtureType
                    << ". Expected 'multiComponent2TMixture', "
                    << "'neReactingMixture' or 'pure2TMixture'. "
                    << "Defaulting mixtureCheck=false." << nl << endl;

                mixtureCheck = false;
            }
        }
        else
        {
            WarningInFunction
                << "thermoType dictionary found, but 'mixture' entry missing. "
                << "Defaulting mixtureCheck=false." << nl << endl;
        }
    }
    else
    {
        WarningInFunction
            << "'thermoType' dictionary missing in thermophysicalProperties. "
            << "Defaulting mixtureCheck=false." << nl << endl;
    }

    // VV relaxation flag
    if (thermoDict_.found("VVrelaxation"))
    {
        useVibVibRelaxation = thermoDict_.get<bool>("VVrelaxation");
        Info<< "VVrelaxation = "
            << (useVibVibRelaxation ? "true" : "false") << nl << endl;
    }
    else
    {
        WarningInFunction
            << "'VVrelaxation' not found in thermophysicalProperties. "
            << "Defaulting to false." << nl << endl;
    }
}

// ── End of thermochemistryCheck.H ────────────────────────────────────────────