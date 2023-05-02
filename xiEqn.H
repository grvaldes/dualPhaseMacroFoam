{
	Info << "Saturation avant: " << " Min(xi) = " << gMin(xi)
	<< " Max(xi) = " << gMax(xi) << endl;

	fvScalarMatrix xiEqn
	(
        	nf*fvm::ddt(xi) + fvc::div(phil)
	);

	xiEqn.solve();

	Info<< "Saturation apres: " << " Min(xi) = " << gMin(xi)
		<< " Max(xi) = " << gMax(xi) << endl;
}
