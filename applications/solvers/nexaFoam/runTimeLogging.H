myFile << "Time is : " << runTime.timeName() << endl;

PtrList<volScalarField>& Q_VT_species = thermo2T.correctVibSource();
PtrList<volScalarField>& Q_VV_species = thermo2T.correctVibVibSource();
PtrList<volScalarField>& Tau_VT_species = thermo2T.correctVTRelaxationTime();

forAll(spec, i)
{
    const volScalarField& Yi = Y_species[i];
    const scalar Wi = mix.W(i);
    const scalar mui = mix.mu(i,TTR.average().value());
    const scalar kappaTRi = mix.kappaTR(i,p.average().value(), TTR.average().value());
    const scalar CpTRi = mix.CpTR(i,p.average().value(), TTR.average().value());

    const scalar volAvg =
        Yi.weightedAverage(mesh.V()).value();

    Info << "Average mass fraction of " << Yi.name() << " is:  " << volAvg << nl;
    myFile << "Average mass fraction of " << Yi.name() << " is:  " << volAvg << nl;

    Info << "Mol weight of " << spec[i] << " is:  " << Wi << nl;
    myFile << "Mol weight of " << spec[i] << " is:  " << Wi << nl;

    Info << "Average viscosity of " << spec[i] << " is:  " << mui << nl;
    myFile << "Average viscosity of " << spec[i] << " is:  " << mui << nl;

    Info << "Average thermal conductivity of " << spec[i] << " is:  " << kappaTRi << nl;
    myFile << "Average thermal conductivity of " << spec[i] << " is:  " << kappaTRi << nl;

    Info << "Average CpTR of " << spec[i] << " is:  " << CpTRi << nl;
    myFile << "Average CpTR of " << spec[i] << " is:  " << CpTRi << nl;

    const volScalarField& Q_VT_i = Q_VT_species[i];

    Info << "Average Q_VT of " << spec[i] << " is: " << Q_VT_i.average().value() << nl;
    myFile << "Average Q_VT of " << spec[i] << " is: " << Q_VT_i.average().value() << nl;

    const volScalarField& Tau_VT_i = Tau_VT_species[i];

    Info << "Average Tau_VT of " << spec[i] << " is: " << Tau_VT_i.average().value() << nl;
    myFile << "Average Tau_VT of " << spec[i] << " is: " << Tau_VT_i.average().value() << nl;

    const volScalarField& Q_VV_i = Q_VV_species[i];

    Info << "Average Q_VV of " << spec[i] << " is: " << Q_VV_i.average().value() << nl;
    myFile << "Average Q_VV of " << spec[i] << " is: " << Q_VV_i.average().value() << nl;

}


myFile << "Average Translational Cv value in internal cells is : \n" << (thermo2T.CvT())().average().value() << endl;
if(vibrationalCheck)
{
    myFile << "Average Rotational Cv value in internal cells is : \n" << (thermo2T.CvR())().average().value() << endl;
}
myFile << "Average Translational-Rotational Cp value in internal cells is : \n" << (thermo2T.CpTR())().average().value() << endl;
if(vibrationalCheck)
{
    myFile << "Average Vibrational Cv value in internal cells is : \n" << (thermo2T.CvVib())().average().value() << endl;
}
myFile << "Average viscosity of value in internal cells is : \n" << thermo2T.wilkeMuAverage() << endl;
Info << "Average viscosity of value in internal cells is : \n" << thermo2T.wilkeMuAverage() << endl;

myFile << "Average thermal conductivity of value in internal cells is : \n" << thermo2T.wilkeKappaAverage() << endl;
Info << "Average thermal conductivity of value in internal cells is : \n" << thermo2T.wilkeKappaAverage() << endl;
if(vibrationalCheck)
{
    myFile << "Average vibrational thermal conductivity value in internal cells is : \n" << (thermo2T.kappaVib())().average().value() << endl;
}

Info << "Average translational-rotational temperature is: " << TTR.average().value() << nl;
myFile << "Average translational-rotational temperature is: " << TTR.average().value() << nl;

if (vibrationalCheck)
{
    Info << "Average vibrational temperature is: " << TVib.average().value() << nl;
    myFile << "Average vibrational temperature is: " << TVib.average().value() << nl;
}

Info << "Average translational-rotational enthalpy is: " << hTR.average().value() << nl;
myFile << "Average translational-rotational enthalpy is: " << hTR.average().value() << nl;

if (vibrationalCheck)
{
    Info << "Average vibrational internal energy is: " << eV.average().value() << nl;
    myFile << "Average vibrational internal energy is: " << eV.average().value() << nl;
}

myFile << "Effective diffusion coefficient of mixture is : \n" << (thermo2T.fickDiffusionCoeff())().average().value() << endl;
Info << "Effective diffusion coefficient of mixture is : \n" << (thermo2T.fickDiffusionCoeff())().average().value() << endl;