{
	Info << "Saturation avant: " << " Min(xi) = " << gMin(xi)
	<< " Max(xi) = " << gMax(xi) << endl;

	Info<< "phig avant = " << average(phig) << nl
        << "phil avant = " << average(phil) << endl;

	fvScalarMatrix xiEqn
	(
        	nf*fvm::ddt(xi) + fvc::div(phil)
	);

	xiEqn.solve();

	Info<< "Saturation apres: " << " Min(xi) = " << gMin(xi)
		<< " Max(xi) = " << gMax(xi) << endl;

	Info<< "phig apres = " << average(phig) << nl
        << "phil apres = " << average(phil) << endl;
}
