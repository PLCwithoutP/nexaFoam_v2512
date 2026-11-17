        // --- Positivity check: stop with a named cell instead of a SIGFPE
        {
            const scalarField& r = rho.primitiveField();
            label badCell = -1;

            forAll(r, celli)
            {
                if (!(r[celli] > 0))        // also catches NaN
                {
                    badCell = celli;
                    break;
                }
            }

            if (returnReduce(badCell >= 0, orOp<bool>()))
            {
                if (badCell >= 0)
                {
                    Pout<< "[positivity] rho = " << r[badCell]
                        << " in cell " << badCell
                        << " at " << mesh.C()[badCell] << endl;
                }

                FatalErrorInFunction
                    << "Non-positive density after the continuity/species "
                    << "update at t = " << runTime.timeName()
                    << ". See the [positivity] line(s) above."
                    << exit(FatalError);
            }
        }
        