if (!mixtureCheck)
{
    solve(fvm::ddt(rho) + fvc::div(phi));
}
else
{
    rhoCont = rho;  // start from current rho
    solve(fvm::ddt(rhoCont) + fvc::div(phi));
}
